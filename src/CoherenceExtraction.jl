export coherence_extraction
export j_out_phase_estimation, initial_state_phase_estimation, fs_pop_phase_estimation
export j_out_compound, coherence_extraction_compound
export fs_pop_compound#, fs_pop_compound_sampled
export j_out_single_setup
#export coherence_extraction_identical, combined_measurement_coherence_extraction_identical
#export j_out4bins, j_out_hom

"""
    coherence_extraction(
        N,
        j_out,
        pops_init,
        pop_fs,
        angles,
        contr_j_tuples=correlated_short_bins_tuples(N, extract_diagonal=false),
        projector_weights=ones(Float64, length(j_out))
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
    from the state. Default value is
    `correlated_short_bins_tuples(N, extract_diagonal=false)`, which is all correlated time
    bins, excluding the correlated populations.
- `projector_weights`: The weights of the projectors in the final state. Default equal
    weights for all projectors.

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
    #ϕ_arr = (0:0.00001:2) * π
    nn_phases = zeros(Float64, N)
    k = 1 # nearest neigbor phase measurements suffice
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
        nn_phases[idx + 1] = -atan(c_imag, c_real)
            # the correction phase(negative sign)
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
        # phases applied through rotated density matrix `ρ_rotated`, not through explicit
        # optional argument `phases` in `explicit_fs_pop`. Should this be changed?

        #pop_fs_real[idx] = explicit_fs_pop(ρ_init, j, angles_real[idx])
        #pop_fs_imag[idx] = explicit_fs_pop(ρ_rotated, j, angles_imag[idx])

    end

    return pop_fs_real, pop_fs_imag
end

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
    n_interference_pairs = floor(Int, N / 2)
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
    pop_contr = contr_pops / N # populations contribution to fidelity

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

    tasks_per_thread = 1
    # customize this as needed. More tasks have more overhead, but better load balancing
    n_settings = length(j_out_all)
    n_threads_compound = min(Base.Threads.nthreads(), 2)
    #chunk_size = max(1, n_settings ÷ (tasks_per_thread * Base.Threads.nthreads()))
    chunk_size = max(1, n_settings ÷ (tasks_per_thread * n_threads_compound))
    data_chunks = Base.Iterators.partition(Array(1:n_settings), chunk_size)
    #for k in eachindex(j_out_all)
    tasks = map(data_chunks) do chunk
        extracted_coherences = 0.0
        Base.Threads.@spawn begin
            for k in chunk
                j_out = j_out_all[k] # the kth final state projector indices
                angles = angles_all[k] # the kth angle settings
                pop_fs = pops_fs_all[k] # the kth final state populations
                contr_j_tuples = contr_j_tuples_all[k]
                # # the two time bins to be extracted from data.

                coh = coherence_extraction(
                    N, j_out, pops_init, pop_fs, angles, contr_j_tuples, projector_weights)
                    # noisy extraction using the correlated and
                    # anti-correlated coincidence counts
                extracted_coherences += coh
            end
            return extracted_coherences
        end
    end
    extracted_coherences_all = sum(fetch.(tasks))
    extracted_coherences_all += pop_contr # add populations contribution to fidelity
    return extracted_coherences_all
end

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

"""
    fs_pop_compound(ρ_init, angles_all=angles_compound(ρ2N(ρ_init)); n_samples=nothing)

Compute all needed final-state populations for an initial state `ρ_init` for a complete set
of measurements in the compound coherence extraction scheme. The `angles_all` argument
contains all beam-splitter angles to be used for the scheme.

Optionally, a number of samples `n_samples` can be specified to compute the final-state
populations with finite sampling statistics.

See also [`fs_pop_phase_estimation`](@ref),
[`explicit_fs_pop`](@ref), [`angles_compound`](@ref).
"""
function fs_pop_compound(ρ_init, angles_all=angles_compound(ρ2N(ρ_init)); n_samples=nothing)

    ρ_init, angles_all, j_out_all, proj_weights = _fs_pop_compound_worker(
        ρ_init, angles_all
    )
    tasks_per_thread = 1
    # customize this as needed. More tasks have more overhead, but better load balancing
    n_settings = length(j_out_all)
    n_threads_compound = min(Base.Threads.nthreads(), 2)
    chunk_size = max(1, n_settings ÷ (tasks_per_thread * n_threads_compound))
    data_chunks = Base.Iterators.partition(Array(1:n_settings), chunk_size)
    tasks = map(data_chunks) do chunk
        Base.Threads.@spawn begin
            pops_out = zeros(Float64, n_settings)
            for k in chunk
                j_out = j_out_all[k]
                angles = angles_all[k]
                pops_out[k] = explicit_fs_pop(
                    ρ_init, j_out, angles, proj_weights; n_samples=n_samples
                ) # Julia 1.3 cannot infere keyword argument name from passed argument name
                # no phase argument here. Not needed, as no phases are being applied iin
                # this step. Density matrix is rotated already in the phase estimation step.
                # In compound scheme, only real parts of some of the coherences are
                # extracted. For imaginary parts of arbitrary coherences, one has to
                # consider phase argument.
            end
            return pops_out
        end
    end
    pops_out_sum = sum(fetch.(tasks))
    return pops_out_sum
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
    #pops_out = Vector{Float64}(undef, length(j_out_all))
    # prepare pops_out array for type stability
    return ρ_init, angles_all, j_out_all, proj_weights#, pops_out

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
