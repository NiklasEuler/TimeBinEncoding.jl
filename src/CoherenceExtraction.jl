export coherence_extraction, compound_coherence_extraction, initial_state_phase_estimation
export angles_kth_neighbor_interference, noisy_angles_symmetric, angles_single_setup
export j_out_single_setup


"""
    coherence_extraction(N, j_out, ρ, angles, noisy_angles=copy(angles), contr_j_idxs = correlated_short_bins_idxs(N); extract_diagonal::Bool=true)

Extract the coherences between correlated time-bin populations.

# Arguments

- `N`: the number of time bins in the initial state `ρ`.
- `j_out`: single `j` index or collection of `j` indices, indicating from which final state projections the coherences get extracted.
- `ρ`: the initial state density matrix from which the coherences are extracted, given in the `|lcmk⟩` basis.
- `angles`: a vector or vectors of scheduled beam-splitter angles. The number of round trips matches `length(angles)`.
- `noisy_angles`: a vector or vectors of actually realized noisy beam-splitter angles used 
    to simulate the measured population. Same shape as `angles`. By default, `noisy_angles`
    is just a copy of `angles`, so no noise is applied.
- `contr_j_idxs`: Vector of `j` indices of states that should be extracted from the state. 
    Default value is `correlated_short_bins_idxs(N)`, which is all correlated time bins.

# Keyword Arguments

- `extract_diagonal`: Bool flag to indicate whether state populations, i.e. diagonal elements of ρ, should be extracted from the coherence. Default is `true`.

See also `compound_coherence_extraction`.
"""
function coherence_extraction(N, j_out, ρ, angles, noisy_angles=copy(angles), contr_j_idxs = correlated_short_bins_idxs(N); extract_diagonal::Bool=true)
    N = convert(Int64, N)::Int64
    j_out = try 
		convert(Vector{Int64}, j_out)::Vector{Int64}
	catch 
		convert(Int64, j_out)::Int64
	end
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    pops = Float64.(diag(ρ))
    j1_arr, j2_arr, weights = explicit_final_state_coherence_map(j_out, angles) # extraction based on assumed ideal angles
    @argcheck weights ≠ []

	pop_j_out_extracted = explicit_final_state_projection_expval(ρ, j_out, noisy_angles) # actual measured value based on noisy angles
	extracted_coherence = []
    extracted_weights = Float64[]
	for idx in eachindex(j1_arr)
		j1 = j1_arr[idx]
		j2 = j2_arr[idx]
		if j1 ∈ contr_j_idxs && j2 ∈ contr_j_idxs && (extract_diagonal || j1 ≠ j2) # check whether both j1 and j2 are correlated time bins
			push!(extracted_coherence,(j1,j2)) # relevant coherence, so indices saved to list of extracted coherences.
            push!(extracted_weights,weights[idx])
		elseif j1 == j2
			pop_j_out_extracted -= pops[j1] * weights[idx] # # non-contributing population. Can be removed exactly as exact value is known.
		else			
			pop_j_out_extracted -= sqrt(pops[j1]*pops[j2]) * abs(weights[idx]) # subtract non-contributing coherence bound.
		end
	end
    @argcheck extracted_weights ≠ []
    n_contr_sched = length(contr_j_idxs)
    n_extracted = length(extracted_weights)
    if (n_extracted != n_contr_sched^2 && extract_diagonal) || (n_extracted != n_contr_sched*(n_contr_sched-1) && !extract_diagonal)
        @warn "Some of the scheduled coherences have a vanishing weight in the given final-state projectors. Please check again and consider adapting the scheduled coherences in `contr_j_idxs`."
    end
    if !all(extracted_weights .≈ extracted_weights[1])
        throw(ArgumentError("The coherences scheduled for extraction have differing weights in the chosen final-state projectors and can thus not straight-forwardly be extracted. Please check input again."))
    else
        norm = extracted_weights[1]
    end
	pop_j_out_extracted /= N*norm # normalization
	return convert(Float64, pop_j_out_extracted)#, extracted_coherence
end


function initial_state_phase_estimation(ρ_init, ϵ_angles=0.0)
    ρ = convert(Matrix{ComplexF64}, copy(ρ_init))::Matrix{ComplexF64}
    ϵ_angles = convert(Float64, ϵ_angles)::Float64
    N = Int64(sqrt(size(ρ)[1]/(n_loops2)))
    
    ϕ_arr = (0:0.00001:2)*π
    nn_phases = zeros(Float64, N)
    k = 1 # nearest neigbour phase measurements suffice
    extract_diagonal = false # dont need the populations
    angles_arr = angles_kth_neighbor_interference(N, k)
    j_out_arr = [[lcmk2j(N+k+1,i,0,i,0),lcmk2j(N+k+1,i+1,1,i+1,1)] for i in 1:k:N-1]
    j_contr_idxs = correlated_short_bins_idxs(N)
    for (idx, j) in enumerate(j_out_arr)
        φ_arr = zeros(Float64, N)
        φ_arr[idx+1] = π/2
        ρ_rotated = phase_on_density_matrix(ρ, φ_arr)
        noisy_angles_real = noisy_angles_symmetric(angles_arr[idx], ϵ_angles)
        noisy_angles_imag = noisy_angles_symmetric(angles_arr[idx], ϵ_angles)
        c_real = coherence_extraction(N, j, ρ, angles_arr[idx], noisy_angles_real, j_contr_idxs[[idx,idx+k]]; extract_diagonal=extract_diagonal)
        c_imag = coherence_extraction(N, j, ρ_rotated, angles_arr[idx], noisy_angles_imag, j_contr_idxs[[idx,idx+k]]; extract_diagonal=extract_diagonal)
        c_contr = c_real .* cos.(-ϕ_arr) .+ c_imag .* sin.(-ϕ_arr)
        nn_phases[idx+1] = ϕ_arr[argmax(c_contr)]
    end
    relative_phases =  mod.(cumsum(nn_phases),2*π)
    ρ_corrected = phase_on_density_matrix(ρ, -1 * relative_phases) # reverse initial state phase profile
    return ρ_corrected, relative_phases
end

function angles_kth_neighbor_interference(N, k)
    N = convert(Int64, N)::Int64
    k = convert(Int64, k)::Int64

    @argcheck N ≥ 1
    @argcheck k ≥ 1
    @argcheck k < N

    angles = Vector{Vector{Vector{Float64}}}(undef, N-k)
    for l in 1:N-k
        angles_measurement = [zeros(Float64, n) for n in N:N+k]
        angles_measurement[1][l] = 0.5 *π
        angles_measurement[end][l+k] = 0.25*π
        angles[l] = angles_measurement
    end
return angles
end

function noisy_angles_symmetric(angles, ϵ_angles)
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    ϵ_angles = convert(Float64, ϵ_angles)::Float64
    noisy_angles = [angles_m .+ ϵ_angles * 2 .* (rand(length(angles_m)) .- 1) for angles_m in angles]
    return noisy_angles
end

function compound_coherence_extraction(ρ, ϵ_angles = 0.0)
    ρ = convert(Matrix{ComplexF64}, copy(ρ))::Matrix{ComplexF64}
    ϵ_angles = convert(Float64, ϵ_angles)::Float64
    N = Int64(sqrt(size(ρ)[1]/(n_loops2)))

    contr_j_idxs_all = correlated_short_bins_idxs(N)
    contr_pops = Float64(sum([ρ[j,j] for j in contr_j_idxs_all]))
    extract_diagonal = false
    extracted_cohereneces = contr_pops/N # populations contribution to fidelity
    for k in 1:N-1
        angles_k = angles_kth_neighbor_interference(N, k)
        j_out_k = [[lcmk2j(N+k+1,i,0,i,0),lcmk2j(N+k+1,i+1,1,i+1,1)] for i in k:1:N-1]
        for (idx, j_out) in enumerate(j_out_k)
            contr_j_idxs = contr_j_idxs_all[[idx,idx+k]]
            angles = angles_k[idx]
            noisy_angles = noisy_angles_symmetric(angles, ϵ_angles)
            coh = coherence_extraction(N, j_out, ρ, angles, noisy_angles, contr_j_idxs; extract_diagonal=extract_diagonal)  # noisy extraction
            extracted_cohereneces += coh
        end
    end
    return extracted_cohereneces
end

function angles_single_setup(N)
    N = convert(Int64, N)::Int64
    @argcheck isinteger(log2(N))
    M = 2*(N-1)
    N_half = N ÷ 2
    pow_half = log2(N_half)
    angles_cascade = [zeros(Float64, n) for n in N:N+M-1]
    angles_cascade[1][1:N_half] .= π/2
    recursive_beam_splitter_array!(N_half, 1 + N_half, 1 + N_half, angles_cascade, "center")
    for i in 2:N_half - 1 # 2 are already included in minimum structure
        angles_cascade[N + i * 2][N + i] = π/4
    end
    return angles_cascade
end

function recursive_beam_splitter_array!(N_bs, n_idx, m_idx, angles, branch)
    @argcheck branch ∈ ["early", "center", "late"]
    
    angles[m_idx][n_idx:n_idx+N_bs-1] .= π/4 # put original bs
    if N_bs == 1 # can happen in the flanks, no further recursion
        if branch == "early"
            angles[m_idx+1][[n_idx]] .= π/2 # left flank transparent couplers
        elseif branch == "late"
            angles[m_idx+1][[n_idx+1]] .= π/2 # right flank transparent couplers
        end
    elseif N_bs == 2 # smallest regular structure in the center, also appears in the flanks. no further recursion 
        angles[m_idx+1][[n_idx,n_idx+2]] .= π/2
        angles[m_idx+1][[n_idx+1]] .= π/4
        angles[m_idx+3][[n_idx+2]] .= π/4
        if branch == "early"
            angles[m_idx+4][[n_idx+1]] .= π/2 # left flank transparent couplers
            angles[m_idx+4][[n_idx+2]] .= π/2 # left flank transparent couplers
        elseif branch == "late"
            angles[m_idx+4][[n_idx+3]] .= π/2 # left flank transparent couplers
            angles[m_idx+4][[n_idx+4]] .= π/2 # left flank transparent couplers
        end
    else
        N_bs_half = Int64(N_bs/2)
        N_bs_quarter = Int64(N_bs/4)
        angles[m_idx+N_bs_half][n_idx:n_idx+N_bs_quarter-1] .= π/2 # left flank transparent couplers
        angles[m_idx+N_bs_half][n_idx+N_bs+N_bs_quarter:n_idx+N_bs+N_bs_half-1] .= π/2 # right flank transparent couplers
        recursive_beam_splitter_array!(N_bs_quarter, n_idx+N_bs_quarter, m_idx+N_bs_half+N_bs_quarter, angles, "early") # left flank beam splitter array
        recursive_beam_splitter_array!(N_bs_quarter, n_idx+N_bs+N_bs_quarter, m_idx+N_bs_half+N_bs_quarter, angles, "late") # right flank beam splitter array
        recursive_beam_splitter_array!(N_bs_half, n_idx+N_bs_half, m_idx+N_bs_half, angles, "center")
    end
end

function j_out_single_setup(N)
    N = convert(Int64, N)::Int64
    @argcheck isinteger(log2(N))
    N_half = Int64(N/2)
	M = 2*(N-1)
	j_short = [lcmk2j(N+M,i,0,i,0) for i in N-1:N-2+N_half]
	j_long = [lcmk2j(N+M,i,1,i,1) for i in N-1+N_half:2*(N-1)]
	j_arr = append!(j_short, j_long)
    return j_arr
end