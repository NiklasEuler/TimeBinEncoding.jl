export coherence_extraction, initial_state_phase_estimation
export angles_kth_neighbor_interference, noisy_angles_symmetric, angles_single_setup
export angles_compound, j_out_compound, compound_coherence_extraction, pops_fs_compound
export j_out_single_setup


"""
    coherence_extraction(
        N, j_out, pops_init, pop_fs, angles,contr_j_idxs = correlated_short_bins_idxs(N);
        extract_diagonal::Bool=true
    )

Extract the coherences between correlated time-bin populations.

# Arguments

- `N`: the number of time bins in the initial state `ρ`.
- `j_out`: single `j` index or collection of `j` indices, indicating from which final state
    projections the coherences get extracted.
- `pops_init`: Vector of the measured initial state population, given in the `|lcmk⟩` basis.
- `pop_fs`: the measured population in the combination of final state projectors given by
    `j_out`
- `angles`: a Vector or Vectors of scheduled beam-splitter angles. The number of round trips
    matches `length(angles)`.
- `contr_j_idxs`: Vector of `j` indices of states that should be extracted from the state.
    Default value is `correlated_short_bins_idxs(N)`, which is all correlated time bins.

# Keyword Arguments

- `extract_diagonal`: Bool flag to indicate whether state populations, i.e., diagonal
    elements of ρ, should be extracted from the coherence. Default is `true`.

See also `compound_coherence_extraction`.
"""
function coherence_extraction(
    N, j_out, pops_init, pop_fs, angles, contr_j_idxs = correlated_short_bins_idxs(N);
    extract_diagonal::Bool=true
)
    N = convert(Int64, N)::Int64
    j_out = try
		convert(Vector{Int64}, j_out)::Vector{Int64}
	catch
		convert(Int64, j_out)::Int64
	end
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    j1_arr, j2_arr, weights = explicit_fs_coherence_map(j_out, angles)
    # extraction based on assumed ideal angles
    @argcheck weights ≠ []

	extracted_coherence = []
    extracted_weights = Float64[]
	for idx in eachindex(j1_arr)
		j1 = j1_arr[idx]
		j2 = j2_arr[idx]
		if j1 in contr_j_idxs && j2 in contr_j_idxs && (extract_diagonal || j1 ≠ j2)
            # check whether both j1 and j2 are correlated time bins
			push!(extracted_coherence, (j1, j2))
            # relevant coherence, so indices saved to list of extracted coherences.
            push!(extracted_weights, weights[idx])
		elseif j1 == j2
			pop_fs -= pops_init[j1] * weights[idx]
            # non-contributing population. Can be removed exactly as exact value is known.
		else
			pop_fs -= sqrt(pops_init[j1]*pops_init[j2]) * abs(weights[idx])
            # subtract non-contributing coherence bound.
		end
	end
    @argcheck extracted_weights ≠ []
    n_contr_sched = length(contr_j_idxs)
    n_extracted = length(extracted_weights)
    if (n_extracted != n_contr_sched^2 && extract_diagonal) ||
         (n_extracted != n_contr_sched * (n_contr_sched - 1) && !extract_diagonal)
        @warn "Some of the scheduled coherences have a vanishing weight in the given "*
        "final-state projectors. Please check again and consider adapting the scheduled "*
        "coherences in `contr_j_idxs`."
    end
    if !all(extracted_weights .≈ extracted_weights[1])
        throw(ArgumentError("The coherences scheduled for extraction have differing "*
        "weights in the chosen final-state projectors and can thus not straight-forwardly "*
        "be extracted. Please check input again."))
    else
        norm = extracted_weights[1]
    end
	pop_fs /= N * norm # normalization
	return convert(Float64, pop_fs)#, extracted_coherence
end

"""
initial_state_phase_estimation(
    ρ_init,
    pops_init=populations(ρ_init),
    angles_real=angles_kth_neighbor_interference(ρ2N(ρ_init), 1),
    angles_imag=angles_kth_neighbor_interference(ρ2N(ρ_init), 1)
)

Compute the relative phases between the time bins of an initial states `ρ_init` with
arbitrary but fixed phases and return the phase-corrected density matrix and measured
relative phases.

Optionally, the measured initial state populations can be given as a parameter, taking
finite sampling statistics into consideration. The default value are the true populations of
the initial state `ρ_init`. Furthermore, the phase estimation process can be modelled to
have noisy beam-splitter angles. The entire noisy angle sets can be given for the real- and
imaginary-part measurements each. The default here are the clean angles.

See also `noisy_angles_symmetric`, `angles_kth_neighbor_interference`,
`phase_on_density_matrix`.
"""
function initial_state_phase_estimation(
    ρ_init,
    pops_init=populations(ρ_init),
    angles_real=angles_kth_neighbor_interference(ρ2N(ρ_init), 1),
    angles_imag=angles_kth_neighbor_interference(ρ2N(ρ_init), 1)
)
    # TODO: Change function to more general noise source, perhaps with noisy angles given
    # directly as a parameter
    ρ = convert(Matrix{ComplexF64}, copy(ρ_init))::Matrix{ComplexF64}
    N = ρ2N(ρ_init) # extract time-bin number from density matrix

    ϕ_arr = (0:0.00001:2) * π
    nn_phases = zeros(Float64, N)
    k = 1 # nearest neigbour phase measurements suffice
    extract_diagonal = false # dont need the populations
    angles_k = angles_kth_neighbor_interference(N, k)
    j_out_arr = [[lcmk2j(N + k + 1, i, 0, i, 0), lcmk2j(N + k + 1, i + 1, 1, i + 1, 1)]
        for i in 1:k:N - 1]
    j_contr_idxs = correlated_short_bins_idxs(N) # all time bin j indices
    for (idx, j) in enumerate(j_out_arr) # j is the index of the final state projectors
        φ_arr = zeros(Float64, N)
        φ_arr[idx + 1] = π / 2
        ρ_rotated = phase_on_density_matrix(ρ, φ_arr)
        pop_fs_real = explicit_fs_projection_expval(ρ_init, j, angles_real[idx])
        pop_fs_imag = explicit_fs_projection_expval(ρ_rotated, j, angles_imag[idx])
        # TODO Implement final state population sampling statistics

        c_real = coherence_extraction(
            N, j, pops_init, pop_fs_real, angles_k[idx], j_contr_idxs[[idx, idx + k]];
            extract_diagonal=extract_diagonal
        )
        c_imag = coherence_extraction(
            N, j, pops_init, pop_fs_imag, angles_k[idx], j_contr_idxs[[idx, idx + k]];
            extract_diagonal=extract_diagonal
        )
        c_contr = c_real .* cos.(-ϕ_arr) .+ c_imag .* sin.(-ϕ_arr)
        nn_phases[idx + 1] = ϕ_arr[argmax(c_contr)]
    end
    relative_phases =  mod.(cumsum(nn_phases), 2 * π)
    ρ_corrected = phase_on_density_matrix(ρ, -1 * relative_phases)
    # reverse initial state phase profile
    return ρ_corrected, relative_phases
end

function ρ2N(ρ)
    return Int64(sqrt(size(ρ)[1] / (N_LOOPS2)))
end

"""
    angles_kth_neighbor_interference(N, k)
    angles_kth_neighbor_interference(N, k, ϵ_angles)

Return an Vector of complete angle settings to exclusively interfere all combination of two
time bins with distance k.

Can also be called with optional argument `ϵ_angles`, which adds uniform random noise to the
angles.

See also `noisy_angles_symmetric`.
"""
function angles_kth_neighbor_interference end

function angles_kth_neighbor_interference(N, k)
    N = convert(Int64, N)::Int64
    k = convert(Int64, k)::Int64

    @argcheck N ≥ 1
    @argcheck k ≥ 1
    @argcheck k < N

    angles = Vector{Vector{Vector{Float64}}}(undef, N - k)
    for i in 1:N - k
        angles[i] = _angles_kth_neighbor(N, k, i)
    end

    return angles
end

function angles_kth_neighbor_interference(N, k, ϵ_angles)
    N = convert(Int64, N)::Int64
    k = convert(Int64, k)::Int64

    @argcheck N ≥ 1
    @argcheck k ≥ 1
    @argcheck k < N

    angles = Vector{Vector{Vector{Float64}}}(undef, N - k)
    for i in 1:N - k
        angles[i] = _angles_kth_neighbor(N, k, i, ϵ_angles)
    end

    return angles
end

function _angles_kth_neighbor(N, k, i)
    angles_k_i = [zeros(Float64, n) for n in N:N + k]
    angles_k_i[1][i] = 0.5  * π # send the early time bin into the long loop
    angles_k_i[end][i+k] = 0.25 * π
    # interfere after k loops / k time bins travelled
    return angles_k_i
end

function _angles_kth_neighbor(N, k, i, ϵ_angles)
    angles_k_i = _angles_kth_neighbor(N, k, i)
    angles_k_i_noisy = noisy_angles_symmetric(angles_k_i, ϵ_angles)
    return angles_k_i_noisy
end

"""
    noisy_angles_symmetric(angles, ϵ_angles)

Given a scheduled beam-splitter angle configuration `angles`, compute a noisy randomly
distributed actual realization of that configuration.

Each noisy angle `̂φ_i` is drawn symmetrically around its planned angle `φ_i`, i.e.,
`̂φ_i ∼ U([φ_i-ϵ_angles, φ_i+ϵ_angles])`.

See also `initial_state_phase_estimation`.
"""
function noisy_angles_symmetric(angles, ϵ_angles)
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    ϵ_angles = convert(Float64, ϵ_angles)::Float64
    noisy_angles = [angles_m .+ ϵ_angles * 2 .* (rand(length(angles_m)) .- 1)
        for angles_m in angles]
    return noisy_angles
end

"""
    compound_coherence_extraction(ρ, ϵ_angles = 0.0)

Extract all correlated time-bin coherences from state `ρ` by surgically interfering all
two-time-bin combinations in a series of different mesh setups.

The process can also be simulated to be noisy by giving the symmetric, uniform random
deviation on the beam-splitter angles `ϵ_angles`.

See also `coherence_extraction`, `noisy_angles_symmetric`.
"""
function compound_coherence_extraction(pops_init, pops_fs_all)
    N = Int64(sqrt(length(pops_init) / N_LOOPS2))::Int64

    contr_j_idxs_all = correlated_short_bins_idxs(N)
    contr_pops = Float64(sum([pops_init[j] for j in contr_j_idxs_all]))
    extract_diagonal = false
    extracted_coherences = contr_pops / N # populations contribution to fidelity
    j_out_all = j_out_compound(N)
    angles_all = angles_compound(N)
    for k in 1:N - 1
        j_out_k = j_out_all[k] # all kth neighbor final state projector indices
        angles_k = angles_all[k] # all kth neighbor angle settings
        pop_fs_k = pops_fs_all[k] # all kth neighbor final state populations
        for (idx, j_out) in enumerate(j_out_k)
            contr_j_idxs = contr_j_idxs_all[[idx, idx + k]]
            # the two time bins to be extracted from data.
            angles = angles_k[idx]
            pop_fs = pop_fs_k[idx]
            coh = coherence_extraction(
                N, j_out, pops_init, pop_fs, angles, contr_j_idxs;
                extract_diagonal = extract_diagonal
            ) # noisy extraction
            extracted_coherences += coh
        end
    end

    return extracted_coherences
end

function j_out_compound(N)
    j_out = [[[lcmk2j(N + k + 1, i, 0, i, 0), lcmk2j(N + k + 1, i + 1, 1, i + 1, 1)]
        for i in k:1:N - 1] for k in 1:N - 1] # pairs of |lSlS⟩ and |l + 1Ll + 1L⟩
    return j_out
end

function angles_compound end

function angles_compound(N)
    angles = [angles_kth_neighbor_interference(N, k) for k in 1:N - 1]
    return angles
end

function angles_compound(N, ϵ_angles)
    angles = [angles_kth_neighbor_interference(N, k, ϵ_angles) for k in 1:N - 1]
    return angles
end

function pops_fs_compound(ρ_init, angles_compound)
    ρ_init = convert(Matrix{ComplexF64}, copy(ρ_init))::Matrix{ComplexF64}
    angles_compound = convert(
        Vector{Vector{Vector{Vector{Float64}}}}, angles_compound
    )::Vector{Vector{Vector{Vector{Float64}}}}
    N = ρ2N(ρ_init)

    j_out_all = j_out_compound(N)
    pops_out = [zeros(Float64, N - k) for k in 1:N - 1]
    for k in 1:N - 1
        j_out_k = j_out_all[k]
        angles_k = angles_compound[k]
        for (idx, j_out) in enumerate(j_out_k)
            angles = angles_k[idx]
            pops_out[k][idx] = explicit_fs_projection_expval(ρ_init, j_out, angles)
        end
    end

    return pops_out
end

"""
    angles_single_setup(N)

Return the angles for a specific beam-splitter setup resulting in full coherence information
in every single readout time bin.

The setup requires `M = 2(N - 1)` round trips and interferes every input state time bin with
every other time bin in each of the non-vanishing output time bins.
"""
function angles_single_setup(N)
    N = convert(Int64, N)::Int64
    @argcheck isinteger(log2(N))
    M = 2 * (N - 1)
    N_half = N ÷ 2
    pow_half = log2(N_half)
    angles_cascade = [zeros(Float64, n) for n in N:N + M - 1]
    angles_cascade[1][1:N_half] .= π / 2
    recursive_beam_splitter_array!(N_half, 1 + N_half, 1 + N_half, angles_cascade, "center")
    for i in 2:N_half - 1 # 2 are already included in minimum structure
        angles_cascade[N + i * 2][N + i] = π / 4
    end

    return angles_cascade
end

function recursive_beam_splitter_array!(N_bs, n_idx, m_idx, angles, branch)
    @argcheck branch in ["early", "center", "late"]

    angles[m_idx][n_idx:n_idx + N_bs - 1] .= π / 4 # put original bs
    if N_bs == 1 # can happen in the flanks, no further recursion
        if branch == "early"
            angles[m_idx + 1][[n_idx]] .= π / 2 # left flank transparent couplers
        elseif branch == "late"
            angles[m_idx + 1][[n_idx + 1]] .= π / 2 # right flank transparent couplers
        end
    elseif N_bs == 2
        # smallest regular structure in the center, also appears in the flanks.
        # no further recursion
        angles[m_idx + 1][[n_idx, n_idx + 2]] .= π / 2
        angles[m_idx + 1][[n_idx + 1]] .= π / 4
        angles[m_idx + 3][[n_idx + 2]] .= π / 4
        if branch == "early"
            angles[m_idx + 4][[n_idx + 1]] .= π / 2 # left flank transparent couplers
            angles[m_idx + 4][[n_idx + 2]] .= π / 2 # left flank transparent couplers
        elseif branch == "late"
            angles[m_idx + 4][[n_idx + 3]] .= π / 2 # left flank transparent couplers
            angles[m_idx + 4][[n_idx + 4]] .= π / 2 # left flank transparent couplers
        end
    else
        N_bs_half = Int64(N_bs / 2)
        N_bs_quarter = Int64(N_bs / 4)
        angles[m_idx + N_bs_half][n_idx:n_idx + N_bs_quarter - 1] .= π / 2
        # left flank transparent couplers
        angles[m_idx + N_bs_half][
            n_idx + N_bs + N_bs_quarter:n_idx + N_bs + N_bs_half - 1] .= π / 2
        # right flank transparent couplers
        recursive_beam_splitter_array!(
            N_bs_quarter, n_idx + N_bs_quarter, m_idx + N_bs_half + N_bs_quarter,
            angles, "early"
        ) # left flank beam splitter array
        recursive_beam_splitter_array!(
            N_bs_quarter, n_idx + N_bs+N_bs_quarter, m_idx + N_bs_half + N_bs_quarter,
            angles, "late"
        ) # right flank beam splitter array
        recursive_beam_splitter_array!(
            N_bs_half, n_idx + N_bs_half, m_idx + N_bs_half, angles, "center"
        )
   end

    return nothing
end

function j_out_single_setup(N)
    N = convert(Int64, N)::Int64
    @argcheck isinteger(log2(N))
    N_half = Int64(N/2)
	M = 2 * (N - 1)
	j_short = [lcmk2j(N + M, i, 0, i, 0) for i in N - 1:N - 2 + N_half]
	j_long = [lcmk2j(N + M, i, 1, i, 1) for i in N - 1 + N_half:2 * (N - 1)]
	j_arr = append!(j_short, j_long)
    return j_arr
end
