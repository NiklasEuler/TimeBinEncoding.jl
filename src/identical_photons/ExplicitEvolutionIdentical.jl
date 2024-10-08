export explicit_fs_projection_identical, explicit_fs_coherence_map_identical
export explicit_fs_pop_identical

function explicit_fs_projection_identical(j_out, angles, phases=ones(Float64, N))
    j_out = convert(Int64, j_out)::Int64 # two-photon bin index in the |l, c , m , k > basis
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    M = length(angles) # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    #if M ≤ 6 # symbolic backend is faster, but too memory intensive for too many iterations
    #    return _explicit_fs_projection_symbolic_backend(N, M, j_out, angles)
    #else
        #return _explicit_fs_projection_mesh_backend(N, M, j_out, angles)
    #end

    #return _explicit_fs_projection_mesh_identical_backend(N, M, j_out, angles, phases)
    return _explicit_fs_projection_mesh_identical_backend_parallel(N, M, j_out, angles, phases)

end

function _explicit_fs_projection_mesh_identical_backend(
    N, M, j_out, angles, phases=ones(Float64, N)
)
    @argcheck abs2.(phases) ≈ ones(Float64, N)
    d_hilbert_space = Int(N_LOOPS * N * (N_LOOPS * N + 1) / 2)
    # full local hilbert space dimension for two photons
    coeff_arr = Vector{ComplexF64}(undef, 0)
    j_idx_arr_contr = Int64[]
    l1, c1, m1, k1, l2, c2, m2, k2 = j_super2lcmk_identical(N + M, j_out)
    # four-photon state indices
#=
    l1_init_min = max(0, l1 - M) # light cone for contributions fron the initial state
    #m1_init_min = max(0, m1 - M)
    l2_init_min = max(0, l2 - M)
    #m2_init_min = max(0, m2 - M)

    l1_init_max = min(N - 1, l1) # light cone for contributions fron the initial state
    m1_init_max = min(N - 1, m1)
    l2_init_max = min(N - 1, l2)
    m2_init_max = min(N - 1, m2)
    # due to symmetry in the two species, could consider only one half of the states
    # would have to consider (j1, j2) and (j2, j1) in the end

    for l1_init in l1_init_min:l1_init_max, m1_init in l1_init:m1_init_max
        for l2_init in l2_init_min:l2_init_max, m2_init in l2_init:m2_init_max =#

    configs_lightcone_arr = _configs_lightcone_identical(N, M, l1, m1, l2, m2)
    for (l1_init, m1_init, l2_init, m2_init) in eachcol(configs_lightcone_arr)
        j_init = lcmk2j_super_identical(N, l1_init, 0, m1_init, 0, l2_init, 0, m2_init, 0)
        single_ket = spzeros(ComplexF64, d_hilbert_space^2)
        phase_idxs = [l1_init, m1_init, l2_init, m2_init] .+ 1 # indexing from 1
        single_ket[j_init] = prod(phases[phase_idxs]) # initial state phase
        single_ket_evolved = mesh_evolution_identical(single_ket, angles)
        coeff = single_ket_evolved[j_out]
        if !isapprox(abs2(coeff), 0.0, atol = WEIGHT_CUTOFF)
            push!(coeff_arr, coeff)
            push!(j_idx_arr_contr, j_init)
        end
    end

    return j_idx_arr_contr, coeff_arr
end

function _explicit_fs_projection_mesh_identical_backend_parallel(
    N, M, j_out, angles, phases=ones(Float64, N)
)
    @argcheck abs2.(phases) ≈ ones(Float64, N)
    #println("Running in parallel!")
    d_hilbert_space = Int(N_LOOPS * N * (N_LOOPS * N + 1) / 2)
    # full local hilbert space dimension for two photons
    #coeff_arr = Vector{ComplexF64}(undef, 0)
    #j_idx_arr_contr = Int64[]
    l1, c1, m1, k1, l2, c2, m2, k2 = j_super2lcmk_identical(N + M, j_out)
    # four-photon state indices
    configs_lightcone_arr = _configs_lightcone_identical(N, M, l1, m1, l2, m2)
    N_comb = size(configs_lightcone_arr, 2)
    tasks_per_thread = 1
    # customize this as needed. More tasks have more overhead, but better load balancing

    chunk_size = max(1, N_comb ÷ (tasks_per_thread * Base.Threads.nthreads()))
    data_chunks = Base.Iterators.partition(Array(1:N_comb), chunk_size)
    # partition your data into chunks that individual tasks will deal with
    # for j in 1:N_comb
    tasks = map(data_chunks) do chunk
        Base.Threads.@spawn begin
            coeff_arr_local = Vector{ComplexF64}(undef, 0)
            j_idx_arr_contr_local = Int64[]
            for j in chunk
                l1_init, m1_init, l2_init, m2_init = configs_lightcone_arr[:, j]
                j_init = lcmk2j_super_identical(
                    N, l1_init, 0, m1_init, 0, l2_init, 0, m2_init, 0
                )
                single_ket = spzeros(ComplexF64, d_hilbert_space^2)
                phase_idxs = [l1_init, m1_init, l2_init, m2_init] .+ 1 # indexing from 1
                single_ket[j_init] = prod(phases[phase_idxs]) # initial state phase
                single_ket_evolved = mesh_evolution_identical(single_ket, angles)
                coeff = single_ket_evolved[j_out]
                if !isapprox(abs2(coeff), 0.0, atol = WEIGHT_CUTOFF)
                    push!(coeff_arr_local, coeff)
                    push!(j_idx_arr_contr_local, j_init)
                end
            end
            return j_idx_arr_contr_local, coeff_arr_local
        end
    end
    results = fetch.(tasks)
	j_idx_arr_contr = cat([tup[1] for tup in results]..., dims = 1)
    coeff_arr = cat([tup[2] for tup in results]..., dims = 1)
    return j_idx_arr_contr, coeff_arr
end

function _configs_lightcone_identical(N, M, l1, m1, l2, m2)
	l1_init_min = max(0, l1 - M) # light cone for contributions fron the initial state
    l2_init_min = max(0, l2 - M)

    l1_init_max = min(N - 1, l1) # light cone for contributions fron the initial state
    m1_init_max = min(N - 1, m1)
    l2_init_max = min(N - 1, l2)
    m2_init_max = min(N - 1, m2)

	N_comb1 = (m1_init_max + 1) * (l1_init_max - l1_init_min + 1) - l1_init_max *
        (l1_init_max + 1) / 2 + (l1_init_min - 1) * l1_init_min / 2
	N_comb2 = (m2_init_max + 1) * (l2_init_max - l2_init_min + 1) - l2_init_max *
        (l2_init_max + 1) / 2 + (l2_init_min - 1) * l2_init_min / 2
	N_comb = Int(N_comb1 * N_comb2)
	configs_lightcone_arr = zeros(Int64, 4, N_comb)
	i = 0
	for l1_init in l1_init_min:l1_init_max, m1_init in l1_init:m1_init_max
        for l2_init in l2_init_min:l2_init_max, m2_init in l2_init:m2_init_max
			i += 1
			configs_lightcone_arr[:, i] = [l1_init, m1_init, l2_init, m2_init]
		end
	end

	return configs_lightcone_arr
end

function explicit_fs_coherence_map_identical end

function explicit_fs_coherence_map_identical(
    j_out::Int64,
    angles,
    projector_weight=1,
    phases::AbstractVector=ones(Float64, length(angles[1]))
)
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    j_idx_arr_contr, coeff_arr = explicit_fs_projection_identical(j_out, angles, phases)
    n_contr = length(j_idx_arr_contr)
    j1_arr = inverse_rle(j_idx_arr_contr, fill(n_contr, n_contr))
    j2_arr = repeat(j_idx_arr_contr, n_contr)
    weights = kron(coeff_arr, conj.(coeff_arr)) * projector_weight
    return j1_arr, j2_arr, weights
end

function explicit_fs_coherence_map_identical(
        j_out_arr::AbstractVector{Int64},
        angles,
        projector_weights=ones(Float64, length(j_out_arr)),
        phases::AbstractVector=ones(Float64, length(angles[1]))
    )
    @argcheck length(j_out_arr) == length(projector_weights)

    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    #M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins
    d_hilbert_space = (N * (2 * N + 1))^2::Int64
        # hilbert space dimension of the initial state
        # simplified expression for N_LOOPS = 2 canceling out geometric series

    j1_arr = Int64[]
    j2_arr = Int64[]
    weights = ComplexF64[]
    weight_vec = SparseVector(d_hilbert_space^2, Int64[], ComplexF64[])

    for (projector_idx, j_out) in enumerate(j_out_arr)
        j_idx_arr_contr, coeff_arr = explicit_fs_projection_identical(j_out, angles, phases)
        for (idx1, j1) in enumerate(j_idx_arr_contr)
            for (idx2, j2) in enumerate(j_idx_arr_contr)
                weight = kron(coeff_arr[idx1], conj(coeff_arr[idx2]))
                j_coh = lm2j(d_hilbert_space, j1 - 1, j2 - 1)
                weight_vec[j_coh] += weight * projector_weights[projector_idx]
            end
        end
    end

    for j_coh in weight_vec.nzind
        j1, j2 = j2lm(d_hilbert_space, j_coh)
        j1 += 1 # correct for base 1 indexing
        j2 += 1 # correct for base 1 indexing
        weight = weight_vec[j_coh]
        if !isapprox(abs2(weight), 0.0, atol = WEIGHT_CUTOFF)
            push!(j1_arr, j1)
            push!(j2_arr, j2)
            push!(weights, weight_vec[j_coh])
        end
   end

    return j1_arr, j2_arr, weights
end

function explicit_fs_pop_identical end

function explicit_fs_pop_identical(
    ρ_init, j_out::Int64, angles, phases::AbstractVector=ones(Float64, length(angles[1]))
)
    projector_weight = 1 # default projector weight
    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map_identical(j_out, angles, projector_weight, phases)
    return expval_calculation(ρ_init, j1_arr, j2_arr, weights)
end

function explicit_fs_pop_identical(
    ρ_init,
    j_out_arr::AbstractVector{Int64},
    angles,
    projector_weights=ones(Float64, length(j_out_arr)),
    phases::AbstractVector=ones(Float64, length(angles[1]))
)
    exp_val = 0.0
    for (j_idx, j_out) in enumerate(j_out_arr)
        exp_val_j = explicit_fs_pop_identical(ρ_init, j_out, angles, phases)
        exp_val_j *= projector_weights[j_idx]
        exp_val += exp_val_j
    end

    return exp_val
end
