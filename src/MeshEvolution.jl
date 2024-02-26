export shift_timebins, beam_splitter_operator, coin_operator, mesh_evolution
export shift_timebins_single_photon
export explicit_ket_evolution_sp, explicit_ket_evolution, explicit_state_evolution
export explicit_final_state_projection_sp, explicit_final_state_projection, explicit_final_state_coherence_map, explicit_add_final_state_projection
export explicit_final_state_projection_expval
export phase_on_density_matrix

global const weight_cutoff = 1e-16

function shift_timebins_single_photon(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    new_vec = Vector{ComplexF64}(undef, length(state_vec)+n_loops)
    new_vec[2] = 0
    new_vec[end-1] = 0
    new_vec[1:2:end-3] = @view state_vec[1:2:end]
    new_vec[4:2:end] = @view state_vec[2:2:end]
    return new_vec
end

function shift_timebins(state_vec::Vector)
    state_vec = convert(Vector{ComplexF64}, state_vec)::Vector{ComplexF64}
    N = Int64(sqrt(length(state_vec)/(n_loops2)))
    new_vec = zeros(ComplexF64, ((N+1)*n_loops)^2)
    for j in eachindex(state_vec)
        l,c,m,k = j2lcmk(N,j)
        shifted_j = lcmk2j(N+1, l+c, c, m+k, k) # adapted system has one more time bin, so we need to put N+1
        new_vec[shifted_j] = state_vec[j]
    end
    return new_vec
end

function beam_splitter_operator(θ)
    θ = convert(Float64,θ)::Float64
    cs = cos(θ)
    sn = im*sin(θ)
    cols = [1,1,2,2]
    rows = [1,2,1,2]
    vals = [cs, sn, sn, cs]
   return sparse(cols,rows,vals)
end

function coin_operator(angles::Vector)
    real_angles = convert(Vector{Float64}, angles)::Vector{Float64}
    matrices = [beam_splitter_operator(θ) for θ in real_angles]
    single_photon_coin_operator = blockdiag(matrices...)
    tensor_coin_operator = kron(single_photon_coin_operator,single_photon_coin_operator)
    return tensor_coin_operator::SparseMatrixCSC{ComplexF64, Int64}
end

"""
    mesh_evolution(ψ_init, angles)

TBW
"""
function mesh_evolution(ψ_init, angles)
    state = convert(Vector{ComplexF64}, ψ_init)::Vector{ComplexF64}
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    for i in eachindex(angles)
        coin_op = coin_operator(angles[i])
        state = coin_op * state
        state = shift_timebins(state)
    end
    return state
end

"""
    symbolic_ket_evolution_sp(M, l)

TBW
"""
function explicit_ket_evolution_sp(l, angles)
    l = convert(Int64, l)::Int64 # initial time bin index
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}} # beam splitter angles
    M = length(angles) # number of roundtrips

    j_idx_arr, trigonometric_history_arr, angle_history_arr = symbolic_ket_evolution_sp(M, l)
    j_idx_arr_contr, coeff_arr = symbolic_2_explicit_worker(angles, j_idx_arr, trigonometric_history_arr, angle_history_arr)
    return j_idx_arr_contr, coeff_arr
end

function explicit_final_state_projection_sp(l, c, angles)
    l = convert(Int64, l)::Int64 # final state time bin index
    c = convert(Int64, c)::Int64 # final state loop index
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}} # beam splitter angles
    M = length(angles) # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs = symbolic_final_state_projection_sp(M, l, c)
    filter_idxs = Int64[]
    for (idx, j) in enumerate(j_idx_arr_fs)
        l,c = j2lc(j)
        if l ≥ N
            push!(filter_idxs, idx) # outside of beam splitter angle range
        end
    end
    deleteat!(j_idx_arr_fs, filter_idxs)
    deleteat!(trigonometric_history_arr_fs, filter_idxs)
    deleteat!(angle_history_arr_fs, filter_idxs)
    j_idx_arr_contr, coeff_arr = symbolic_2_explicit_worker(angles, j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs)
    return j_idx_arr_contr, coeff_arr
end

function symbolic_2_explicit_worker(angles, j_idx_arr, trigonometric_history_arr, angle_history_arr)
    M = length(angles) # number of roundtrips
    coeff_arr = Vector{ComplexF64}(undef, 0)
    j_idx_arr_contr = Int64[]
    tri_vals = [hcat(cos.(ang), sin.(ang)) for ang in angles]

    for (i, j) in enumerate(j_idx_arr)
        coeff = Complex(0)
        for k in 1:size(angle_history_arr[i])[1]
            tri_string = trigonometric_history_arr[i][k,:]
            phase_factor = im^(sum(tri_string))
            angle_string = angle_history_arr[i][k,:] .+ 1
            tri_string .+= 1 # shift all values to 1 and 2, making them their respective indices of tri_vals matrices
            coeff += prod([tri_vals[m][angle_string[m], tri_string[m]] for m in 1:M]) * phase_factor
        end
        if !isapprox(abs2(coeff), 0.0, atol=weight_cutoff)
            push!(coeff_arr, coeff)
            push!(j_idx_arr_contr, j)
        end
    end
    return j_idx_arr_contr, coeff_arr
end

function explicit_ket_evolution(j_init, angles)
    j_init = convert(Int64, j_init)::Int64 # two-photon bin index in the |l,c,m,k> basis

    M = length(angles) # number of roundtrips
    N = length(angles[1]) # initial number of time bins
    n_bins = N + M # maximum number of bins after evolution

    l_init, c_init, m_init, k_init = j2lcmk(N, j_init)
    @argcheck c_init == 0
	@argcheck k_init == 0
    j_idx_arr_l, coeff_arr_l = explicit_ket_evolution_sp(l_init, angles)
    j_idx_arr_m, coeff_arr_m = explicit_ket_evolution_sp(m_init, angles)
    j_idx_arr_contr, coeff_arr = sp_2_two_photon(n_bins, j_idx_arr_l, j_idx_arr_m, coeff_arr_l, coeff_arr_m)

    return j_idx_arr_contr, coeff_arr
end

function explicit_final_state_projection(j_out, angles)
    j_out = convert(Int64, j_out)::Int64 # two-photon bin index in the |l,c,m,k> basis

    M = length(angles) # number of roundtrips
    N = length(angles[1]) # initial number of time bins
    n_bins = N # maximum number of bins before evolution ≕ N

    l_out, c_out, m_out, k_out = j2lcmk(N+M, j_out)

    j_idx_arr_l, coeff_arr_l = explicit_final_state_projection_sp(l_out, c_out, angles)
    j_idx_arr_m, coeff_arr_m = explicit_final_state_projection_sp(m_out, k_out, angles)
    j_idx_arr_contr, coeff_arr = sp_2_two_photon(n_bins, j_idx_arr_l, j_idx_arr_m, coeff_arr_l, coeff_arr_m)
    
    return j_idx_arr_contr, coeff_arr
end

function sp_2_two_photon(n_bins, j_idx_arr_l, j_idx_arr_m, coeff_arr_l, coeff_arr_m)
    coeff_arr = Vector{ComplexF64}(undef, 0)
    j_idx_arr_contr = Int64[]

    for (idxl, jl) in enumerate(j_idx_arr_l)
        l,c = j2lc(jl)
        for (idxm, jm) in enumerate(j_idx_arr_m)
            m,k = j2lc(jm)
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

    Ψ_out = zeros(ComplexF64, n_loops2*(N+M)^2)
    for (j_init, coeff) in enumerate(state)
		if coeff == 0
			continue
		end
        j_idx_arr_contr, coeff_arr = explicit_ket_evolution(j_init, angles)
        Ψ_out[j_idx_arr_contr] .+= coeff .* coeff_arr
    end
    return Ψ_out
end


function explicit_final_state_coherence_map(j_out::Int64, angles)
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    j_idx_arr_contr, coeff_arr = explicit_final_state_projection(j_out, angles)

    n_contr = length(j_idx_arr_contr)
    j1_arr = inverse_rle(j_idx_arr_contr,fill(n_contr,n_contr))
    j2_arr = repeat(j_idx_arr_contr,n_contr)
    weights = kron(coeff_arr,conj.(coeff_arr))
    return j1_arr, j2_arr, weights
end

function explicit_final_state_coherence_map(j_out_arr::Vector{Int64}, angles)
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    #M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins
    d_hilbert_space = n_loops2*N^2::Int64 # hilbert space dimension of the initial state

    j1_arr = Int64[]
    j2_arr = Int64[]
    weights = ComplexF64[]
    weight_vec = SparseVector(d_hilbert_space^2,Int64[],ComplexF64[])
    for j_out in j_out_arr
        j_idx_arr_contr, coeff_arr = explicit_final_state_projection(j_out, angles)
        for (idx1, j1) in enumerate(j_idx_arr_contr), (idx2, j2) in enumerate(j_idx_arr_contr)
            weight = kron(coeff_arr[idx1],conj(coeff_arr[idx2]))
            j_coh = lm2j(d_hilbert_space, j1, j2)
            weight_vec[j_coh] += weight
        end
    end
    for j_coh in weight_vec.nzind
        j1, j2 = j2lm(d_hilbert_space, j_coh)
        weight = weight_vec[j_coh]
        if !isapprox(abs2(weight), 0.0, atol = weight_cutoff)  
            push!(j1_arr, j1)
            push!(j2_arr, j2)
            push!(weights, weight_vec[j_coh])
        end
    end
    return j1_arr, j2_arr, weights
end

function explicit_final_state_projection_expval(ρ_init, j_out::Int64, angles)
    j1_arr, j2_arr, weights = explicit_final_state_coherence_map(j_out, angles)
    exp_val = 0
    for i in eachindex(j1_arr)
        j1 = j1_arr[i]
        j2 = j2_arr[i]
        weight = weights[i]
        exp_val += ρ_init[j1,j2]*weight
    end
    return convert(Float64, real(exp_val))
end

function explicit_final_state_projection_expval(ρ_init, j_out_arr::Vector{Int64}, angles)
    exp_val = 0.0
    for j_out in j_out_arr
        exp_val += explicit_final_state_projection_expval(ρ_init, j_out, angles)
    end
    return exp_val
end

function phase_on_density_matrix(ρ, φ_arr)
    ρ_rot = convert(Matrix{ComplexF64}, copy(ρ))::Matrix{ComplexF64}
    φ_arr = convert(Vector{Float64}, φ_arr)::Vector{Float64}
    N = Int64(sqrt(size(ρ)[1]/(n_loops2)))

    @argcheck length(φ_arr) == N
    
    tb_idxs = [lcmk2j(N,i,0,i,0) for i in 0:N-1]
    for (idx1, j1) in enumerate(tb_idxs), (idx2, j2) in enumerate(tb_idxs)
        ρ_rot[j1,j2] *= cis(φ_arr[idx1]-φ_arr[idx2])
    end
    return ρ_rot
end