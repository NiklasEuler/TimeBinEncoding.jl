export coherence_extraction
export j_out_phase_estimation, initial_state_phase_estimation, fs_pop_phase_estimation
export j_out_compound, coherence_extraction_compound
export fs_pop_compound, fs_pop_compound_sampled
export j_out_single_setup
#export coherence_extraction_identical, combined_measurement_coherence_extraction_identical
#export j_out4bins, j_out_hom

"""
    coherence_extraction(
        N, j_out, pops_init, pop_fs, angles, contr_j_tuples=correlated_short_bins_tuples(N);
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
- `contr_j_tuples`: Vector of `(j1, j2)` tuples of state indices that should be extracted
    from the state. Default value is `correlated_short_bins_tuples(N)`, which is all
    correlated time bins.

# Keyword Arguments

- `extract_diagonal`::Bool: Bool flag to indicate whether state populations, i.e., diagonal
    elements of ρ, should be extracted from the coherence. Default is `true`.

See also [`coherence_extraction_compound`](@ref).
"""
function coherence_extraction(
    N,
    j_out,
    pops_init,
    pop_fs,
    angles,
    contr_j_tuples = correlated_short_bins_tuples(N, extract_diagonal=false),
    projector_weights=ones(Float64, length(j_out))
)
    N = convert(Int64, N)::Int64
    j_out = try
		convert(Vector{Int64}, j_out)::Vector{Int64}
	catch
		convert(Int64, j_out)::Int64
	end
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}

    if length(j_out) == 1
        projector_weights = projector_weights[1]
    end

    projector_weights = try
		convert(Vector{Float64}, projector_weights)::Vector{Float64}
	catch
		convert(Float64, projector_weights)::Float64
	end
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}

    j1_arr, j2_arr, weights = explicit_fs_coherence_map(j_out, angles, projector_weights)
    # extraction based on assumed ideal angles
    @argcheck weights ≠ []

	extracted_coherence = []
    extracted_weights = Float64[]

    for idx in eachindex(j1_arr)
		j1 = j1_arr[idx]
		j2 = j2_arr[idx]
		if (j1, j2) in contr_j_tuples
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

    n_contr_sched = length(contr_j_tuples)
    n_extracted = length(extracted_weights)

    if n_extracted != n_contr_sched
        @warn "Some of the scheduled coherences have a vanishing weight in the given "*
        "final-state projectors. Please check again and consider adapting the scheduled "*
        "coherences in `contr_j_tuples`."
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
    initial_state_phase_estimation(pops_init, pops_fs_real, pops_fs_imag)

Compute the relative phases between the time bins of an initial states `ρ_init` with
arbitrary but fixed phases and return the phase-corrected density matrix and measured
relative phases.

The function takes the measured initial-state and final-state populations as an argument,
taking finite sampling statistics and noisy beam-splitter angles into consideration.

See also [`noisy_angles_symmetric`](@ref), [`angles_phase_estimation`](@ref),
[`phase_on_density_matrix`](@ref).
"""
function initial_state_phase_estimation(pops_init, pops_fs_real, pops_fs_imag)
    N = Int64(sqrt(length(pops_init) / N_LOOPS2))
    # extract time-bin number from population number
    ϕ_arr = (0:0.00001:2) * π
    nn_phases = zeros(Float64, N)
    k = 1 # nearest neigbor phase measurements suffice
    extract_diagonal = false # dont need the populations
    #angles_k = angles_kth_neighbor_interference(N, k)
    angles = angles_phase_estimation(N) # intended angles (no noise)
    j_out_arr = j_out_phase_estimation(N)
    j_contr_idxs = correlated_short_bins_idxs(N) # all time bin j indices

    for (idx, j) in enumerate(j_out_arr) # j is the index of the final state projectors
        j_contr1 = j_contr_idxs[idx]
        j_contr2 = j_contr_idxs[idx + k]
        j_contr = [(j_contr1, j_contr2), (j_contr2, j_contr1)]
        c_real = coherence_extraction(
            N, j, pops_init, pops_fs_real[idx], angles, j_contr)
        c_imag = coherence_extraction(
            N, j, pops_init, pops_fs_imag[idx], angles, j_contr)
        #c_real = coherence_extraction(
        #    N, j, pops_init, pops_fs_real[idx], angles_k[idx], j_contr_idxs[[idx, idx + k]];
        #    extract_diagonal=extract_diagonal
        #)
        #c_imag = coherence_extraction(
        #    N, j, pops_init, pops_fs_imag[idx], angles_k[idx], j_contr_idxs[[idx, idx + k]];
        #    extract_diagonal=extract_diagonal
        #)
        #println("idx: ", idx, " lcmk: ", [j2lcmk(N + 2, j[1]), j2lcmk(N + 2, j[2])])
        #println("c_real: ", c_real)
        #println("c_imag: ", c_imag)
        c_contr = c_real .* cos.(-ϕ_arr) .+ c_imag .* sin.(-ϕ_arr)
        nn_phases[idx + 1] = ϕ_arr[argmax(c_contr)]
    end

    relative_phases =  mod.(cumsum(nn_phases), 2 * π)
    return relative_phases
end

"""
    j_out_phase_estimation(N)

Return the vector of all necessary `j` indices for a phase-estimation measurement scheme of
an initial state with `N` time bins.

The returned vector contains vectors with two `j` entries each, corresponding to the two
kets `|i,S,i,S⟩` and `|i+1,L,i+1,L⟩`, for i in 1:N-1.

See also [`j_out_compound`](@ref), [`j_out_single_setup`](@ref), [angles]
[`initial_state_phase_estimation`](@ref), [`fs_pop_phase_estimation`](@ref).
"""
function j_out_phase_estimation(N)
    k = 1 # nearest neigbour phase measurements suffice
    j_out_arr = [[lcmk2j(N + k + 1, i, 0, i, 0), lcmk2j(N + k + 1, i + 1, 1, i + 1, 1)]
    # + k + 1 because two round trips are needed for nearest-neighbor interference
        for i in 1:k:N - k]
    return j_out_arr
end

"""
    function fs_pop_phase_estimation(
        ρ_init,
        angles_real=angles_phase_estimation(ρ_init),
        angles_imag=angles_phase_estimation(ρ_init),
    )

Compute the final state populations measured for the iniital state `ρ_init` in an optionally
noisy setup, given through the two sets of optional mesaurement angles `angles_real` and
`angles_imag`.

Their default values are the noiseless angles for nearest-neighbor time-bin interference.

See also [`fs_pop_compound`](@ref), [`initial_state_phase_estimation`](@ref),
[`angles_phase_estimation`](@ref).
"""
function fs_pop_phase_estimation(
        ρ_init,
        angles_real=angles_phase_estimation(ρ_init),
        angles_imag=angles_phase_estimation(ρ_init),
)
    N = ρ2N(ρ_init) # extract N from ρ_init

    j_out_arr = j_out_phase_estimation(N)
    pop_fs_real = zeros(Float64, length(j_out_arr))
    pop_fs_imag = zero(pop_fs_real)

    φ_arr = (0:N - 1) .* (π / 2)
    ρ_rotated = phase_on_density_matrix(ρ_init, φ_arr)

    for (idx, j) in enumerate(j_out_arr) # j is the index of the final state projectors

        #φ_arr = zeros(Float64, N)
        #φ_arr[idx + 1] = π / 2 # apply π / 2 phase shift to swap real and imaginary parts

        #ρ_rotated = phase_on_density_matrix(ρ_init, φ_arr)


        pop_fs_real[idx] = explicit_fs_pop(ρ_init, j, angles_real)
        pop_fs_imag[idx] = explicit_fs_pop(ρ_rotated, j, angles_imag)

        #pop_fs_real[idx] = explicit_fs_pop(ρ_init, j, angles_real[idx])
        #pop_fs_imag[idx] = explicit_fs_pop(ρ_rotated, j, angles_imag[idx])

    end

    return pop_fs_real, pop_fs_imag
end


#= """
    coherence_extraction_compound(pops_init, pops_fs_all)

Extract all correlated time-bin coherences from the intial- and final-state populations
`pops_init` and `pops_fs_all` by surgically interfering all two-time-bin combinations in a
series of different mesh setups.

See also [`coherence_extraction`](@ref), [`noisy_angles_symmetric`](@ref).
"""
function coherence_extraction_compound(pops_init, pops_fs_all)
    N = Int64(sqrt(length(pops_init) / N_LOOPS2))::Int64
    extract_diagonal = false


    contr_j_idxs_all = correlated_short_bins_idxs(N)
    contr_pops = Float64(sum([pops_init[j] for j in contr_j_idxs_all]))
    extracted_coherences = contr_pops / N # populations contribution to fidelity

    j_out_all = j_out_compound(N)
    angles_all = angles_compound(N)
    projector_weights = [1, 1, -1, -1] # correlated outcomes minues anti-correlated outcomes

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
                N, j_out, pops_init, pop_fs, angles, contr_j_idxs, projector_weights;
                extract_diagonal = extract_diagonal
            ) # noisy extraction using the correlated and anti-correlated coincidence counts
            extracted_coherences += coh
        end
    end

    return extracted_coherences
end
 =#


"""
    proj_weights_compound(N)

Compute the projector weights for the compound-system readout scheme. Correlated outcomes
have weight 1, wherease anti-correlated outcomes have weight -1.

# Arguments
- `N`: The number of particles in the system.

# Returns
- `Vector{Int}`: Array of projector weights.

"""
function proj_weights_compound(N)
    n_interference_pairs = ceil(Int, N / 2)
    projector_weights = [1, 1, -1, -1] # correlated outcomes minues anti-correlated outcomes
    projector_weights_flex =
        repeat(projector_weights, n_interference_pairs) # projector weights
    return projector_weights_flex
end

"""
    coherence_extraction_compound(pops_init, pops_fs_all)

Extract all correlated time-bin coherences from the intial- and final-state populations
`pops_init` and `pops_fs_all` by surgically interfering all two-time-bin combinations in a
series of different mesh setups.

See also [`coherence_extraction`](@ref), [`noisy_angles_symmetric`](@ref).
"""
function coherence_extraction_compound(pops_init, pops_fs_all)
    N = Int64(sqrt(length(pops_init) / N_LOOPS2))::Int64

    contr_j_idxs = correlated_short_bins_idxs(N)
    contr_pops = Float64(sum([pops_init[j] for j in contr_j_idxs]))
    extracted_coherences = contr_pops / N # populations contribution to fidelity

    pairings = graph_coloring(N)

    contr_j_tuples_all = Vector{Vector{Tuple{Int64, Int64}}}(undef, length(pairings))
    # construct vector of all relevant coherences to be extracted from current iteration
    for (idx, pairs) in enumerate(pairings)
        contr_j_tuples_all[idx] = []
        for pair in pairs
            j1 = contr_j_idxs[pair[1] + 1]
            j2 = contr_j_idxs[pair[2] + 1]
            push!(contr_j_tuples_all[idx], (j1, j2), (j2, j1))
        end
    end


    j_out_all = j_out_compound(N)
    angles_all = angles_compound(N)
    projector_weights = proj_weights_compound(N)


    for k in eachindex(j_out_all)
        j_out = j_out_all[k] # the kth final state projector indices
        angles = angles_all[k] # the kth angle settings
        pop_fs = pops_fs_all[k] # the kth final state populations
        contr_j_tuples = contr_j_tuples_all[k]
        # # the two time bins to be extracted from data.

        coh = coherence_extraction(
            N, j_out, pops_init, pop_fs, angles, contr_j_tuples, projector_weights)
            # noisy extraction using the correlated and anti-correlated coincidence counts
        extracted_coherences += coh
    end

    return extracted_coherences
end


 #= """
    j_out_compound(N)

Return all the final-state projector indices for the compound coherence-extraction scheme.

At the lowest layer, always four indices are bundeled together in a Vector, corresponding to
|i,S,i,S⟩ and |i + 1,L,i + 1,L⟩ and the two mixtures |i,S,i+1,L⟩ and |i+1,L,i,S⟩,
respectively. One layer up, all of these Vector corresponding to the intereference of time
bins with distance `k` in the initial state are grouped in a Vector. Finally, all such
grouping for values `k` from 1 to N - 1 are collected in the outermost layer.

See also [`angles_compound`](@ref), [`coherence_extraction_compound`](@ref),
[`j_out_single_setup`](@ref).
"""
function j_out_compound(N)
    j_out = [
        [
            [
                lcmk2j(N + k + 1, i, 0, i, 0),
                lcmk2j(N + k + 1, i + 1, 1, i + 1, 1),
                lcmk2j(N + k + 1, i, 0, i + 1, 1),
                lcmk2j(N + k + 1, i + 1, 1, i, 0)
            ]
            for i in k:1:N - 1
        ]
        for k in 1:N - 1
    ] # pairs of |i,S,i,S⟩ and |i + 1,L,i + 1,L⟩ and mixtures |i,S,i+1,L⟩ and |i+1,L,i,S⟩
    return j_out
end =#

 """
    j_out_compound(N)

Return all the final-state projector indices for the compound coherence-extraction scheme.

# Returns
- `Vector{Vector{Int}}`: All final-state projector indices for each beam-splitter setup.

See also [`angles_compound`](@ref), [`coherence_extraction_compound`](@ref),
[`j_out_single_setup`](@ref).
"""
function j_out_compound(N)
    j_out_arr = Vector{Int}[]
    pairings = graph_coloring(N)
    proj_bin_idxs = [[pairing[j][2] for j in eachindex(pairing)] for pairing in pairings]
    for (idx, pair_arr) in enumerate(pairings)
        j_out = Int[]
        Δpairs = [pair[2] - pair[1] for pair in pair_arr]
        Δmax = max(Δpairs...) # maximum distance between time-bin pairs
        M = Δmax + 1 # number of round trips
        proj_bins = proj_bin_idxs[idx]
        for j in proj_bins
            push!(j_out,
                lcmk2j(N + M, j, 0, j, 0), # correlated time bins
                lcmk2j(N + M, j + 1, 1, j + 1, 1), # correlated time bins
                lcmk2j(N + M, j, 0, j + 1, 1), # anti-correlated time bins
                lcmk2j(N + M, j + 1, 1, j, 0) # anti-correlated time bins
            )
        end
        # pairs of |i,S,i,S⟩ and |i + 1,L,i + 1,L⟩ and mixtures |i,S,i+1,L⟩ and |i+1,L,i,S⟩
        push!(j_out_arr, j_out)
    end
    return j_out_arr
end


#= """
    fs_pop_compound(ρ_init, angles_all=angles_compound(ρ2N(ρ_init)))

Compute all needed final-state populations for an initial state `ρ_init` for a complete set
of measurements in the compound coherence extraction scheme. The `angles_all` argument
contains all beam-splitter angles to be used for the scheme.

See also [`fs_pop_phase_estimation`](@ref), [`explicit_fs_pop`](@ref),
[`angles_compound`](@ref).
"""
function fs_pop_compound(ρ_init, angles_all=angles_compound(ρ2N(ρ_init)))
    # Add projector weights to explicit_fs_pop
    ρ_init = convert(Matrix{ComplexF64}, copy(ρ_init))::Matrix{ComplexF64}
    angles_all = convert(
        Vector{Vector{Vector{Vector{Float64}}}}, angles_all
    )::Vector{Vector{Vector{Vector{Float64}}}}
    N = ρ2N(ρ_init)

    j_out_all = j_out_compound(N)
    pops_out = [zeros(Float64, N - k) for k in 1:N - 1]

    for k in 1:N - 1
        j_out_k = j_out_all[k]
        angles_k = angles_all[k]
        for (idx, j_out) in enumerate(j_out_k)
            angles = angles_k[idx]
            pops_out[k][idx] = explicit_fs_pop(ρ_init, j_out, angles)
        end
    end

    return pops_out
end =#

"""
    fs_pop_compound(ρ_init, angles_all=angles_compound(ρ2N(ρ_init)))

Compute all needed final-state populations for an initial state `ρ_init` for a complete set
of measurements in the compound coherence extraction scheme. The `angles_all` argument
contains all beam-splitter angles to be used for the scheme.

See also [`fs_pop_compound_sampled`](@ref), [`fs_pop_phase_estimation`](@ref),
[`explicit_fs_pop`](@ref), [`angles_compound`](@ref).
"""
function fs_pop_compound(ρ_init, angles_all=angles_compound(ρ2N(ρ_init)))

    ρ_init, angles_all, j_out_all, proj_weights, pops_out = _fs_pop_compound_worker(
        ρ_init, angles_all
    )
    for k in eachindex(j_out_all)
        j_out = j_out_all[k]
        angles = angles_all[k]
        pops_out[k] = explicit_fs_pop(ρ_init, j_out, angles, proj_weights)
    end
    return pops_out
end

"""
    fs_pop_compound_sampled(ρ_init, n_samples, angles_all=angles_compound(ρ2N(ρ_init)))

Compute the sampled final-state populations for an initial state `ρ_init` for a complete set
of measurements in the compound coherence extraction scheme, using `n_sample` sta. The
`angles_all` argument contains all beam-splitter angles to be used for the scheme.


See also [`fs_pop_compound`](@ref), [`fs_pop_phase_estimation`](@ref),
[`explicit_fs_pop`](@ref), [`angles_compound`](@ref).
"""
function fs_pop_compound_sampled(
    ρ_init, n_samples, angles_all=angles_compound(ρ2N(ρ_init))
)
    ρ_init, angles_all, j_out_all, proj_weights, pops_out = _fs_pop_compound_worker(
        ρ_init, angles_all
    )
    for k in eachindex(j_out_all)
        j_out = j_out_all[k]
        angles = angles_all[k]
        pops_out[k] = explicit_fs_pop_sampled(
            ρ_init, j_out, angles, n_samples, proj_weights
        )
    end
    return pops_out
end


function _fs_pop_compound_worker(ρ_init, angles_all)
    # Add projector weights to explicit_fs_pop
    ρ_init = convert(Matrix{ComplexF64}, copy(ρ_init))::Matrix{ComplexF64}
    angles_all = convert(
        Vector{Vector{Vector{Float64}}}, angles_all
    )::Vector{Vector{Vector{Float64}}}
    N = ρ2N(ρ_init)

    j_out_all = j_out_compound(N)
    proj_weights = proj_weights_compound(N) # projector weights
    pops_out = Vector{Float64}(undef, length(j_out_all))

    return ρ_init, angles_all, j_out_all, proj_weights, pops_out

end

"""
    j_out_single_setup(N)

Return all the final-state projector indices for the single-setup coherence-extraction
scheme.

See also ['j_out_compound'](@ref), ['j_out_phase_estimation'](@ref).
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


#= function coherence_extraction_identical(
    N,
    j_out,
    pops_init,
    pop_fs,
    angles,
    contr_j_tuples=correlated_short_bins_tuples_identical(N),
    projector_weights=ones(Float64, length(j_out));
    extract_diagonal::Bool=true
)
    N = convert(Int64, N)::Int64
    j_out = try
		convert(Vector{Int64}, j_out)::Vector{Int64}
	catch
		convert(Int64, j_out)::Int64
	end

    if length(j_out) == 1
        projector_weights = projector_weights[1]
    end

    projector_weights = try
		convert(Vector{Float64}, projector_weights)::Vector{Float64}
	catch
		convert(Float64, projector_weights)::Float64
	end
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}

    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map_identical(j_out, angles, projector_weights)
        # extraction based on assumed ideal angles
    @argcheck weights ≠ []

	extracted_coherence = []
    extracted_weights = Float64[]
    for idx in eachindex(j1_arr)
		j1 = j1_arr[idx]
		j2 = j2_arr[idx]
		if (j1,j2) in contr_j_tuples && (extract_diagonal || j1 ≠ j2)
            # check whether both j1 and j2 are correlated time bins
			push!(extracted_coherence, (j1, j2))
            # relevant coherence, so indices saved to list of extracted coherences.
            push!(extracted_weights, weights[idx])
		elseif j1 == j2
			pop_fs -= pops_init[j1] * weights[idx]
            # non-contributing population. Can be removed exactly as exact value is known.
		else
			pop_fs -= sqrt(pops_init[j1] * pops_init[j2]) * abs(weights[idx])
            if sqrt(pops_init[j1] * pops_init[j2]) * abs(weights[idx]) != 0

                #println(
                    #"Subtracting: (",j1, " ", j2, ") ", sqrt(pops_init[j1]*pops_init[j2]) *
                    #abs(weights[idx])
                #)
            end
            # subtract non-contributing coherence bound.
		end
	end

    @argcheck extracted_weights ≠ []

    n_contr_sched = length(contr_j_tuples)
    n_extracted = length(extracted_weights)

    if (n_extracted != n_contr_sched )
        @warn "Some of the scheduled coherences have a vanishing weight in the given "*
        "final-state projectors. Please check again and consider adapting the scheduled "*
        "coherences in `contr_j_tuples`."
    end

    if !all(extracted_weights .≈ extracted_weights[1])
        @warn "The coherences scheduled for extraction have differing "*
        "weights in the chosen final-state projectors and can thus not straight-forwardly "*
        "be extracted. Please check input again."
        norm = 1 #minimum(extracted_weights) #maximum(extracted_weights)
        # solution for now, but not ideal.
        # Better solution would compute sensible normalization.
    else
        norm = extracted_weights[1]
    end


	pop_fs /= norm # normalization of the extracted coherences

	return convert(Float64, pop_fs)#, extracted_coherence
end

function j_out4bins(N_post, bin_1, bin_2)
    j_SSSS = lcmk2j_super_identical(N_post, bin_1, 0, bin_1, 0, bin_1, 0, bin_1, 0)
    j_LLLL = lcmk2j_super_identical(N_post, bin_2, 1, bin_2, 1, bin_2, 1, bin_2, 1)
    j_SSLL = lcmk2j_super_identical(N_post, bin_1, 0, bin_1, 0, bin_2, 1, bin_2, 1)
    j_LLSS = lcmk2j_super_identical(N_post, bin_2, 1, bin_2, 1, bin_1, 0, bin_1, 0)
    j_SLSL = lcmk2j_super_identical(N_post, bin_1, 0, bin_2, 1, bin_1, 0, bin_2, 1)

    j_out = [j_SSSS, j_LLLL, j_SSLL, j_LLSS, j_SLSL]
    return j_out
end

function j_out_hom(N_post, bin_1, bin_2)
    j_out_lm_SS = lcmk2j_super_identical(N_post, bin_1, 0, bin_1, 0, bin_1, 0, bin_1, 0)
    j_out_lm_LL = lcmk2j_super_identical(N_post, bin_2, 1, bin_2, 1, bin_2, 1, bin_2, 1)
    j_out_lm_hom = lcmk2j_super_identical(N_post, bin_1, 0, bin_2, 1, bin_1, 0, bin_2, 1,)

    j_out = [j_out_lm_SS, j_out_lm_LL, j_out_lm_hom]
    return j_out
end

function combined_measurement_coherence_extraction_identical(
    N,
    combined_weights,
    pops_init,
    pop_fs_combined,
    contr_j_tuples=correlated_short_bins_tuples_identical(N);
    extract_diagonal::Bool=true
)
    d_full_hilbert_space = lcmk2j_super_identical(N, N-1, 1, N-1, 1, N-1, 1, N-1, 1)
	#extracted_coherence = []
    pop_fs = copy(pop_fs_combined)
    extracted_weights = Float64[]
    for j in combined_weights.nzind
		j1, j2 = j2lm(d_full_hilbert_space, j)
		j1 += 1
		j2 += 1
		if (j1,j2) in contr_j_tuples && (extract_diagonal || j1 ≠ j2)
            # check whether both j1 and j2 are correlated time bins
			#push!(extracted_coherence, (j1, j2))
            # relevant coherence, so indices saved to list of extracted coherences.
            push!(extracted_weights, combined_weights[j])
		elseif j1 == j2
			pop_fs -= pops_init[j1] * combined_weights[j]
            # non-contributing population. Can be removed exactly as exact value is known.
		else
			pop_fs -= sqrt(pops_init[j1] * pops_init[j2]) * abs(combined_weights[j])
            # subtract non-contributing coherence bound.
		end
	end

    @argcheck extracted_weights ≠ []

    n_contr_sched = length(contr_j_tuples)
    n_extracted = length(extracted_weights)

    if (n_extracted != n_contr_sched )
        @warn "Some of the scheduled coherences have a vanishing weight in the given "*
        "final-state projectors. Please check again and consider adapting the scheduled "*
        "coherences in `contr_j_tuples`."
    end

    if !all(extracted_weights .≈ extracted_weights[1])
        @warn "The coherences scheduled for extraction have differing "*
        "weights in the chosen final-state projectors and can thus not straight-forwardly "*
        "be extracted. Please check input again."
        norm = 1 #minimum(extracted_weights) #maximum(extracted_weights)
        # solution for now, but not ideal.
        # Better solution would compute sensible normalization.
    else
        norm = extracted_weights[1]
    end


	pop_fs /= norm # normalization of the extracted coherences

	return convert(Float64, pop_fs)#, extracted_coherence
end =#
