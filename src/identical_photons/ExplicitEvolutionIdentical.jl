export explicit_fs_projection_sp_identical, explicit_fs_coherence_map_identical
export explicit_fs_pop_identical


"""
    explicit_fs_projection_sp_identical(j_out, angles, phases)

Compute the backwards propagation of two indistinguishable photons in the state with basis
index `j_out` through the loop system with beam splitter `angles` and initial-state
`phases`.

# Arguments
- `j_out`: Two-photon bin index in the |l, c, m, k⟩ basis.
- `angles`: Vector of beam splitter angles for each roundtrip.
- `phases`: Vector of initial-state phase factors for each time bin.

# Returns
- `j_idx_arr_contr`: Vector of basis indices that contribute to the final state of `j_out`.
- `coeff_arr`: Vector of coefficients for each contributing basis index.

See also [`explicit_fs_coherence_map_identical`](@ref), [`explicit_fs_pop_identical`](@ref).

"""
function explicit_fs_projection_sp_identical(
    j_out, angles, phases=ones(Float64, length(angles[1]))
)
    j_out = convert(Int64, j_out)::Int64 # two-photon bin index in the |l, c , m , k > basis
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    M = length(angles) # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    return _explicit_fs_projection_mesh_identical_backend(
        N, M, j_out, angles, phases
    )

end

"""
    _kron_mem_arr(N, M)

Preallocate memory for the Kronecker products of the coin operators. Return a vector with
an preallocated sparse array for each roundtrip.

See also [`_explicit_fs_projection_mesh_identical_backend`](@ref).
"""
function _kron_mem_arr(N, M)
    dims = [(n * (2 * n + 1)) ^ 2 for n in N:N + M]
    mem = [spzeros(ComplexF64, dim, dim) for dim in dims]
    return mem
end

function _explicit_fs_projection_mesh_identical_backend(
    N, M, j_out, angles, phases=ones(Float64, N)
)
    @argcheck abs2.(phases) ≈ ones(Float64, N)
    #println("Running in parallel!")
    d_hilbert_space_init = Int(N * (2 * N + 1))
    # full local hilbert space dimension for two photons before evolution / after backwards
    # propagation
    d_hilbert_space_final = Int((N + M) * (2 * (N + M) + 1))
    # full local hilbert space dimension for two photons after evolution / before backwards
    # propagation
    coeff_arr = ComplexF64[]
    j_idx_arr_contr = Int64[]

    kron_mem = _kron_mem_arr(N, M) # preallocate memory for kron products
    coin_op = [adjoint(coin_operator_identical(angles[i], kron_mem[i])) for i in M:-1:1]
    # adjoint coin operators for backwards propagation
    shift_op = [adjoint(shift_timebins_operator_identical(i)) for i in N+M-1:-1:N]
    # adjoint shift operators for backwards propagation
    state = spzeros(ComplexF64, d_hilbert_space_final^2)
    state[j_out] = 1
    # final state with singular contribution from j_out
    phases_inv = conj(phases)
    for i in eachindex(coin_op)
        state = shift_op[i] * state
        state = coin_op[i] * state
    end
    for j in 1:d_hilbert_space_init^2
        l1, c1, m1, k1, l2, c2, m2, k2 = j_super2lcmk_identical(N, j)

        coeff = state[j]
        if !isapprox(abs2(coeff), 0.0, atol = WEIGHT_CUTOFF)
            # find nonzero contributions
            phase_idxs = [l1, m1, l2, m2] .+ 1 # indexing from 1
            coeff *= prod(phases_inv[phase_idxs]) # inverse initial state phase
            push!(coeff_arr, coeff)
            push!(j_idx_arr_contr, j)
        end
    end
    return j_idx_arr_contr, coeff_arr
end

"""
    explicit_fs_coherence_map_identical(
        j_out::Real, angles, projector_weight=1, phases
    )
    explicit_fs_coherence_map_identical(
        j_out_arr::AbstractVector{<:Real}, angles, projector_weights, phases
    )

Compute the list of all coherences that contribute to the full final state of two
indistinguishable photon species after propagation through the loop system with beam-
splitter `angles` and initial-state `phases`. Also accepts a `projector_weight` argument to
scale the contributuing coherences. A single final state `j_out` or a vector of final states
`j_out_arr` can be provided. In the latter case, the sum of all contributions, modified by
the corresponding `projector weights` is returned.

# Arguments
- `j_out`: four-photon bin index in the |l1, c1, m1, k1, l2, c2, m2, k2⟩ basis.
    Can be a single integer or a vector of integers. In the second case, the contributions
    are summed up with weights according to `projector_weights`.
- `angles`: Vector of beam splitter angles for each roundtrip.
- `projector_weight`: Weight for the projector contribution. For an integer `j_out`, the
    projector contribution is simply scaled by this factor. For a vector `j_out_arr`, the
    respective contributions are scaled by the corresponding element of `projector_weights`.
- `phases`: Vector of initial-state phase factors for each time bin.

# Returns
- `j1_arr`: Vector of first basis indices of the contributing coherences in the
    |l1, c1, m1, k1, l2, c2, m2, k2⟩ basis.
- `j2_arr`: Vector of second basis indices of the contributing coherences in the
    |l1, c1, m1, k1, l2, c2, m2, k2⟩ basis.
- `weights`: Vector of weights for each contributing coherence.

See also [`explicit_fs_projection_sp_identical`](@ref), [`explicit_fs_pop_identical`](@ref).

"""
function explicit_fs_coherence_map_identical end

function explicit_fs_coherence_map_identical(
    j_out::Real,
    angles,
    projector_weight=1,
    phases::AbstractVector=ones(Float64, length(angles[1]))
)

    projector_weight = convert(Float64, projector_weight)::Float64
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    j_idx_arr_contr, coeff_arr = explicit_fs_projection_sp_identical(j_out, angles, phases)
    n_contr = length(j_idx_arr_contr)
    j1_arr = inverse_rle(j_idx_arr_contr, fill(n_contr, n_contr))::Vector{Int64}
    j2_arr = repeat(j_idx_arr_contr, n_contr)::Vector{Int64}
    weights = kron(coeff_arr, conj.(coeff_arr)) * projector_weight
    return j1_arr, j2_arr, weights::Vector{ComplexF64}
end

function explicit_fs_coherence_map_identical(
        j_out_arr::AbstractVector{<:Real},
        angles,
        projector_weights=ones(Float64, length(j_out_arr)),
        phases::AbstractVector=ones(Float64, length(angles[1]))
    )
    @argcheck length(j_out_arr) == length(projector_weights)

    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    projector_weights = convert(Vector{Float64}, projector_weights)::Vector{Float64}

    N = length(angles[1]) # initial number of time bins
    d_hilbert_space = (N * (2 * N + 1))^2::Int64
        # hilbert space dimension of the initial state
        # simplified expression for N_LOOPS = 2 canceling out geometric series

    j1_arr = Int64[]
    j2_arr = Int64[]
    weights = ComplexF64[]
    weight_vec = SparseVector(d_hilbert_space^2, Int64[], ComplexF64[])

    for (projector_idx, j_out) in enumerate(j_out_arr)
        j_idx_arr_contr, coeff_arr =
            explicit_fs_projection_sp_identical(j_out, angles, phases)
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


"""
    explicit_fs_pop_identical(
        ρ_init, j_out::Real, angles, phases
    )
    explicit_fs_pop_identical(
        ρ_init, j_out_arr::AbstractVector{<:Real}, angles, projector_weights, phases
    )

Compute the expectation value of the populations of `j_out`/`j_out_arr` in the full four-
photon state `ρ_init` after evolution through the loop system with beam-splitter `angles`
and initial-state `phases`. Either, a single final state `j_out` or a vector of final states
can be provided. In the latter case, the sum of all populations, modified by the
corresponding `projector weights` is returned.

# Arguments
- `ρ_init`: Initial density matrix of the four-photon state.
- `j_out`: four-photon bin index in the |l1, c1, m1, k1, l2, c2, m2, k2⟩ basis.
    Can be a single integer or a vector of integers. In the latter case, the populations
    are summed up with weights according to `projector_weights`.
- `angles`: Vector of beam splitter angles for each roundtrip.
- `projector_weights`: Weights for the projector contribution. Only for a vector `j_out_arr`.
    The respective contributions are scaled by the corresponding element of
    `projector_weights`.
- `phases`: Vector of initial-state phase factors for each time bin.

# Returns
- `exp_val`: Expectation value of the populations of `j_out`/`j_out_arr` in the full four-
    photon state `ρ_init`.

See also [`explicit_fs_projection_sp_identical`](@ref),
[`explicit_fs_coherence_map_identical`](@ref).

"""
function explicit_fs_pop_identical end

function explicit_fs_pop_identical(
    ρ_init, j_out::Real, angles, phases::AbstractVector=ones(Float64, length(angles[1]))
)
    projector_weight = 1 # default projector weight
    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map_identical(j_out, angles, projector_weight, phases)
    return expval_calculation(ρ_init, j1_arr, j2_arr, weights)
end

function explicit_fs_pop_identical(
    ρ_init,
    j_out_arr::AbstractVector{<:Real},
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
