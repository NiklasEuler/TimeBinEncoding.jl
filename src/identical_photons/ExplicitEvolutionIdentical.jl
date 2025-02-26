export explicit_fs_projection_identical, explicit_fs_coherence_map_identical
export explicit_fs_pop_identical



function explicit_fs_projection_identical(
    j_out, angles, phases=ones(Float64, length(angles[1]))
)
    j_out = convert(Int64, j_out)::Int64 # two-photon bin index in the |l, c , m , k > basis
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    M = length(angles) # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    return _explicit_fs_projection_mesh_identical_backend_parallel(
        N, M, j_out, angles, phases
    )

end

function _kron_mem_arr(N, M)
    dims = [(n * (2 * n + 1)) ^ 2 for n in N:N + M]
    mem = [spzeros(ComplexF64, dim, dim) for dim in dims]
    return mem
end

function _explicit_fs_projection_mesh_identical_backend_parallel(
    N, M, j_out, angles, phases=ones(Float64, N)
)
    @argcheck abs2.(phases) ≈ ones(Float64, N)
    #println("Running in parallel!")
    d_hilbert_space_init = Int(N * (2 * N + 1))
    d_hilbert_space_final = Int((N + M) * (2 * (N + M) + 1))
    # full local hilbert space dimension for two photons
    coeff_arr = ComplexF64[]
    j_idx_arr_contr = Int64[]


    kron_mem = _kron_mem_arr(N, M) # preallocate memory for kron products
    coin_op = [adjoint(coin_operator_identical(angles[i], kron_mem[i])) for i in M:-1:1]
    shift_op = [adjoint(shift_timebins_operator_identical(i)) for i in N+M-1:-1:N]
    state = spzeros(ComplexF64, d_hilbert_space_final^2)
    state[j_out] = 1
    phases_inv = conj(phases)
    for i in eachindex(coin_op)
        state = shift_op[i] * state
        state = coin_op[i] * state
    end
    for j in 1:d_hilbert_space_init^2
        l1, c1, m1, k1, l2, c2, m2, k2 = j_super2lcmk_identical(N, j)

        coeff = state[j]
        if !isapprox(abs2(coeff), 0.0, atol = WEIGHT_CUTOFF)
            phase_idxs = [l1, m1, l2, m2] .+ 1 # indexing from 1
            coeff *= prod(phases_inv[phase_idxs]) # initial state phase
            push!(coeff_arr, coeff)
            push!(j_idx_arr_contr, j)
        end
    end
    return j_idx_arr_contr, coeff_arr
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
    j1_arr = inverse_rle(j_idx_arr_contr, fill(n_contr, n_contr))::Vector{Int64}
    j2_arr = repeat(j_idx_arr_contr, n_contr)::Vector{Int64}
    weights = kron(coeff_arr, conj.(coeff_arr)) * projector_weight
    return j1_arr, j2_arr, weights::Vector{ComplexF64}
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
