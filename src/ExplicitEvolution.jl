export explicit_ket_evolution_sp, explicit_ket_evolution, explicit_state_evolution
export explicit_fs_projection_sp, explicit_fs_projection,
    explicit_fs_coherence_map, explicit_add_fs_projection
export explicit_fs_projection_expval


"""
    explicit_ket_evolution_sp(l, angles)

Evolve the single-photon state `|l, 0⟩' according to `angles`, using a symbolic backend.

# Output
- `j_idx_arr_contr`: vector of `j` indices with nonzero contribution in the state after
    evolution.
- `coeff_arr`: vector of corresponding coefficients to the `j` indices given in
    `j_idx_arr_contr`.

See also  [`explicit_ket_evolution`](@ref), [`explicit_fs_projection`](@ref),
[`explicit_fs_projection_sp`](@ref).
"""
function explicit_ket_evolution_sp(l, angles)
    l = convert(Int64, l)::Int64 # initial time bin index
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    # beam splitter angles
    M = length(angles) # number of roundtrips

    j_idx_arr, trigonometric_history_arr, angle_history_arr =
        symbolic_ket_evolution_sp(M, l)
    j_idx_arr_contr, coeff_arr = _symbolic_2_explicit_worker(
            angles, j_idx_arr, trigonometric_history_arr, angle_history_arr
        )
    return j_idx_arr_contr, coeff_arr
end

"""
    explicit_fs_projection_sp(l, c, angles)

Computes the indicies and weights of all initial states that contribute to the single-photon
state `|lc⟩` after evolution by `angles`.

# Output
- `j_idx_arr_contr`: vector of `j` indices with nonzero contribution to `|l, c ⟩` after
    evolution.
- `coeff_arr`: vector of corresponding coefficients to the `j` indices given in
    `j_idx_arr_contr`.

See also [`explicit_fs_projection`](@ref), [`explicit_ket_evolution_sp`](@ref),
[`explicit_fs_coherence_map`](@ref).
"""
function explicit_fs_projection_sp(l, c, angles)
    l = convert(Int64, l)::Int64 # final state time bin index
    c = convert(Int64, c)::Int64 # final state loop index
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    # beam splitter angles
    M = length(angles) # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs =
        symbolic_fs_projection_sp(M, l, c)
    filter_idxs = Int64[]
    for (idx, j) in enumerate(j_idx_arr_fs)
        l, c  = j2lc(j)
        if l ≥ N
            push!(filter_idxs, idx) # outside of beam splitter angle range
        end
    end
    deleteat!(j_idx_arr_fs, filter_idxs)
    deleteat!(trigonometric_history_arr_fs, filter_idxs)
    deleteat!(angle_history_arr_fs, filter_idxs)
    j_idx_arr_contr, coeff_arr = _symbolic_2_explicit_worker(
        angles, j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs
    )
    return j_idx_arr_contr, coeff_arr
end

function _symbolic_2_explicit_worker(
    angles, j_idx_arr, trigonometric_history_arr, angle_history_arr
)
    M = length(angles) # number of roundtrips
    coeff_arr = Vector{ComplexF64}(undef, 0)
    j_idx_arr_contr = Int64[]
    tri_vals = [hcat(cos.(ang), sin.(ang)) for ang in angles]

    for (i, j) in enumerate(j_idx_arr)
        coeff = Complex(0)
        for k in 1:size(angle_history_arr[i])[1]
            tri_string = trigonometric_history_arr[i][k, :]
            phase_factor = im^(sum(tri_string))
            angle_string = angle_history_arr[i][k, :] .+ 1
            # shift all time bins by one, such that the first (l = 0) time bin aligns
            # with the index 1
            tri_string .+= 1
            # shift all values to 1 and 2, making them their respective indices of
            # tri_vals matrices
            coeff += (prod([tri_vals[m][angle_string[m], tri_string[m]] for m in 1:M])
            * phase_factor)
        end
        if !isapprox(abs2(coeff), 0.0, atol = WEIGHT_CUTOFF)
            push!(coeff_arr, coeff)
            push!(j_idx_arr_contr, j)
            # only keep idxs of contributing coefficients, i.e. coeff ≠ 0
        end
   end

    return j_idx_arr_contr, coeff_arr
end

"""
    explicit_ket_evolution(j_init, angles)

Evolve the two-photon state related to `j` according to `angles`, using a symbolic backend.

# Output
- `j_idx_arr_contr`: vector of `j` indices with nonzero contribution in the state after
    evolution.
- `coeff_arr`: vector of corresponding coefficients to the `j` indices given in
    `j_idx_arr_contr`.

See also [`explicit_ket_evolution_sp`](@ref), [`mesh_evolution`](@ref),
[`explicit_fs_projection`](@ref).
"""
function explicit_ket_evolution(j_init, angles)
    j_init = convert(Int64, j_init)::Int64
    # two-photon bin index in the |l, c , m , k > basis

    M = length(angles) # number of roundtrips
    N = length(angles[1]) # initial number of time bins
    n_bins = N + M # maximum number of bins after evolution

    l_init, c_init, m_init, k_init = j2lcmk(N, j_init)
    @argcheck c_init == 0
	@argcheck k_init == 0
    j_idx_arr_l, coeff_arr_l = explicit_ket_evolution_sp(l_init, angles)
    j_idx_arr_m, coeff_arr_m = explicit_ket_evolution_sp(m_init, angles)
    j_idx_arr_contr, coeff_arr =
        _sp_2_two_photon(n_bins, j_idx_arr_l, j_idx_arr_m, coeff_arr_l, coeff_arr_m)

    return j_idx_arr_contr, coeff_arr
end

function explicit_fs_projection(j_out, angles)
    j_out = convert(Int64, j_out)::Int64 # two-photon bin index in the |l, c , m , k > basis
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    M = length(angles) # number of roundtrips
    N = length(angles[1]) # initial number of time bins
    if M ≤ 6 # symbolic backend is faster, but too memory intensive for too many iterations
        return _explicit_fs_projection_symbolic_backend(N, M, j_out, angles)
    else
        return _explicit_fs_projection_mesh(N, M, j_out, angles)
    end
end

function _explicit_fs_projection_symbolic_backend(N, M, j_out, angles)
    l_out, c_out, m_out, k_out = j2lcmk(N + M, j_out)

    j_idx_arr_l, coeff_arr_l = explicit_fs_projection_sp(l_out, c_out, angles)
    j_idx_arr_m, coeff_arr_m = explicit_fs_projection_sp(m_out, k_out, angles)
    j_idx_arr_contr, coeff_arr =
        _sp_2_two_photon(N, j_idx_arr_l, j_idx_arr_m, coeff_arr_l, coeff_arr_m)
    return j_idx_arr_contr, coeff_arr
end

function _explicit_fs_projection_mesh(N, M, j_out, angles)
    coeff_arr = Vector{ComplexF64}(undef, 0)
    j_idx_arr_contr = Int64[]
    l_out, c_out, m_out, k_out = j2lcmk(N + M, j_out)
    if c_out == 0
        l_init_max = min(N - 1, l_out)
        # light cone for contributions fron the initial state
    else
        l_init_max = min(N - 1, l_out - 1)
        # |n,S> cannot contribute to |n,L>, at most |n-1,S> can
    end
    if k_out == 0
        m_init_max = min(N - 1, m_out)
        # light cone for contributions fron the initial state
    else
        m_init_max = min(N - 1, m_out - 1)
        # |n,S> cannot contribute to |n,L>, at most |n-1,S> can
    end
    for l_init in 0:l_init_max, m_init in 0:m_init_max
        j_init = lcmk2j(N, l_init, 0, m_init, 0)
        #single_ket = zeros(ComplexF64, N_LOOPS2 * N^2)
        single_ket = spzeros(ComplexF64, N_LOOPS2 * N^2)
        single_ket[j_init] = 1.0
        single_ket_evolved = mesh_evolution(single_ket, angles)
        coeff = single_ket_evolved[j_out]
        if !isapprox(abs2(coeff), 0.0, atol = WEIGHT_CUTOFF)
            push!(coeff_arr, coeff)
            push!(j_idx_arr_contr, j_init)
        end
   end

    return j_idx_arr_contr, coeff_arr
end

function _sp_2_two_photon(n_bins, j_idx_arr_l, j_idx_arr_m, coeff_arr_l, coeff_arr_m)
    coeff_arr = Vector{ComplexF64}(undef, 0)
    j_idx_arr_contr = Int64[]
    for (idxl, jl) in enumerate(j_idx_arr_l)
        l, c  = j2lc(jl)
        for (idxm, jm) in enumerate(j_idx_arr_m)
            m, k  = j2lc(jm)
            j = lcmk2j(n_bins, l, c, m, k)
            coeff = coeff_arr_l[idxl] * coeff_arr_m[idxm]
            if coeff ≠ 0
                push!(coeff_arr, coeff)
                push!(j_idx_arr_contr, j)
            end
        end
   end

    return j_idx_arr_contr, coeff_arr
end

function explicit_state_evolution(Ψ_init, angles)
    state = convert(Vector{ComplexF64}, Ψ_init)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins
    Ψ_out = zeros(ComplexF64, N_LOOPS2 * (N + M)^2)
    for (j_init, coeff) in enumerate(state)
		if coeff == 0
			continue
		end
        j_idx_arr_contr, coeff_arr = explicit_ket_evolution(j_init, angles)
        Ψ_out[j_idx_arr_contr] .+= coeff .* coeff_arr
   end

    return Ψ_out
end


function explicit_fs_coherence_map(j_out::Int64, angles)
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    j_idx_arr_contr, coeff_arr = explicit_fs_projection(j_out, angles)
    n_contr = length(j_idx_arr_contr)
    j1_arr = inverse_rle(j_idx_arr_contr, fill(n_contr, n_contr))
    j2_arr = repeat(j_idx_arr_contr, n_contr)
    weights = kron(coeff_arr, conj.(coeff_arr))
    return j1_arr, j2_arr, weights
end

function explicit_fs_coherence_map(j_out_arr::Vector{Int64}, angles)
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    #M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins
    d_hilbert_space = N_LOOPS2 * N^2::Int64 # hilbert space dimension of the initial state
    j1_arr = Int64[]
    j2_arr = Int64[]
    weights = ComplexF64[]
    weight_vec = SparseVector(d_hilbert_space^2, Int64[], ComplexF64[])
    for j_out in j_out_arr
        j_idx_arr_contr, coeff_arr = explicit_fs_projection(j_out, angles)
        for (idx1, j1) in enumerate(j_idx_arr_contr)
            for (idx2, j2) in enumerate(j_idx_arr_contr)
                weight = kron(coeff_arr[idx1], conj(coeff_arr[idx2]))
                j_coh = lm2j(d_hilbert_space, j1, j2)
                weight_vec[j_coh] += weight
            end
        end
    end
    for j_coh in weight_vec.nzind
        j1, j2 = j2lm(d_hilbert_space, j_coh)
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
    explicit_fs_projection_expval(ρ_init, j_out::Int64, angles)
    explicit_fs_projection_expval(ρ_init, j_out_arr::Vector{Int64}, angles)


TBW
"""
function explicit_fs_projection_expval end

function explicit_fs_projection_expval(ρ_init, j_out::Int64, angles)
    j1_arr, j2_arr, weights = explicit_fs_coherence_map(j_out, angles)
    return expval_calculation(ρ_init, j1_arr, j2_arr, weights)
end

function explicit_fs_projection_expval(ρ_init, j_out_arr::Vector{Int64}, angles)
    exp_val = 0.0
    for j_out in j_out_arr
        exp_val += explicit_fs_projection_expval(ρ_init, j_out, angles)
   end

    return exp_val
end

function expval_calculation(ρ_init, j1_arr, j2_arr, weights)
    exp_val = 0.0
    for i in eachindex(j1_arr)
        j1 = j1_arr[i]
        j2 = j2_arr[i]
        weight = weights[i]
        exp_val += ρ_init[j1, j2] * weight
   end

    return convert(Float64, real(exp_val))
end
