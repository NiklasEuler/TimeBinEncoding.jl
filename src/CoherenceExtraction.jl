export coherence_extraction, initial_state_phase_estimation
export angles_kth_neighbor_interference, noisy_angles_symmetric, angles_single_setup
export angles_compound, j_out_compound, compound_coherence_extraction, pops_fs_compound
export j_out_single_setup


"""
    coherence_extraction(
        N, j_out, pops_init, pop_fs, angles, contr_j_idxs = correlated_short_bins_idxs(N);
        extract_diagonal::Bool=true
    )

Extract the coherences between correlated time-bin populations.

# Arguments

- `N`: the number of time bins in the initial state `ρ`.
- `j_out`: single `j` index or collection of `j` indices, indicating from which final state
    projections the coherences get extracted.
- `pops_init`: Vector of the measured initial state population, given in the `|lcmk⟩` basis.
- `pop_fs`: Vector of measured final state population in the combination of final state
    projectors given by `j_out`
- `angles`: a Vector or Vectors of scheduled beam-splitter angles. The number of round trips
    matches `length(angles)`.
- `contr_j_idxs`: Vector of `j` indices of states that should be extracted from the state.
    Default value is `correlated_short_bins_idxs(N)`, which is all correlated time bins.

# Keyword Arguments

- `extract_diagonal`::Bool: Bool flag to indicate whether state populations, i.e., diagonal
    elements of ρ, should be extracted from the coherence. Default is `true`.

See also [`compound_coherence_extraction`](@ref).
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

See also [`noisy_angles_symmetric`](@ref), [`angles_kth_neighbor_interference`](@ref),
[`phase_on_density_matrix`](@ref).
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
        φ_arr[idx + 1] = π / 2 # apply π / 2 phase shift to swap real and imaginary parts
        ρ_rotated = phase_on_density_matrix(ρ, φ_arr)

        pop_fs_real = explicit_fs_projection_expval(ρ_init, j, angles_real[idx])
        pop_fs_imag = explicit_fs_projection_expval(ρ_rotated, j, angles_imag[idx])

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

"""
    compound_coherence_extraction(pops_init, pops_fs_all)

Extract all correlated time-bin coherences from the intial- and final-state populations
`pops_init` and `pops_fs_all` by surgically interfering all two-time-bin combinations in a
series of different mesh setups.

See also [`coherence_extraction`](@ref), [`noisy_angles_symmetric`](@ref).
"""
function compound_coherence_extraction(pops_init, pops_fs_all)
    N = Int64(sqrt(length(pops_init) / N_LOOPS2))::Int64
    extract_diagonal = false


    contr_j_idxs_all = correlated_short_bins_idxs(N)
    contr_pops = Float64(sum([pops_init[j] for j in contr_j_idxs_all]))
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

"""
    j_out_compound(N)

Return all the final-state projector indices for the compound coherence-extraction scheme.

At the lowest layer, always two indices are bundeled together in a Vector, corresponding to
|i,S,i,S⟩ and |i + 1,L,i + 1,L⟩, respectively. One layer up, all of these Vector
corresponding to the intereference of time bins with distance `k` in the initial state are
grouped in a Vector. Finally, all such grouping for values `k` from 1 to N - 1 are collected
in the outermost layer.

See also [`angles_compound`](@ref), [`compound_coherence_extraction`](@ref),
[`j_out_single_setup`](@ref).
"""
function j_out_compound(N)
    j_out = [[[lcmk2j(N + k + 1, i, 0, i, 0), lcmk2j(N + k + 1, i + 1, 1, i + 1, 1)]
        for i in k:1:N - 1] for k in 1:N - 1] # pairs of |i,S,i,S⟩ and |i + 1,L,i + 1,L⟩
    return j_out
end

"""
    pops_fs_compound(ρ_init, angles_compound)

Compute all needed final-state populations for an initial state `ρ_init` for a complete set
of measurements in the compound coherence extraction scheme. The `angles_compound` argument
contains all beam-splitter angles to be used for the scheme.

See also [`explicit_fs_projection_expval`](@ref), [`angles_compound`](@ref).
"""
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
    j_out_single_setup(N)

Return all the final-state projector indices for the single-setup coherence-extraction
scheme.

See also ['j_out_single_setup'](@ref).
"""
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
