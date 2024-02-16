export shift_timebins, beam_splitter_operator, coin_operator, mesh_evolution
export correlated_timebin_state, insert_initial_state
export shift_timebins_single_photon
export explicit_ket_evolution_sp

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

function correlated_timebin_state(wf_coeffs::Vector)#todo: normalize
    wf_coeffs = convert(Vector{ComplexF64}, wf_coeffs)::Vector{ComplexF64}
    N = length(wf_coeffs)
    coeffs = normalize(wf_coeffs)
    time_bin_state_vec = zeros(ComplexF64, N^2)
    for i in 0:N-1
        j = lm2j(N,i,i)
        time_bin_state_vec[j] = coeffs[i+1]
    end
    return time_bin_state_vec
end


function insert_initial_state(time_bin_state_vec::Vector)
    N = Int64(sqrt(length(time_bin_state_vec)))
    full_state_vec = zeros(ComplexF64, length(time_bin_state_vec)*n_loops2)
    for l in 0:N-1, m in 0:N-1
        j_time_bin = lm2j(N, l, m)
        j_full = lcmk2j(N, l, 0, m, 0) #insertion into the short-short ket
        full_state_vec[j_full] = time_bin_state_vec[j_time_bin]
    end
    return full_state_vec
end

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
function explicit_ket_evolution_sp(M, l, angles)
    M = convert(Int64, M)::Int64 # number of roundtrips
    l = convert(Int64, l)::Int64 # initial time bin index
    #@argcheck M > 0
    @argcheck l ≥ 0
    j_idx_arr, trigonometric_history_arr, angle_history_arr = symbolic_ket_evolution_sp(M, l)
    coeff_arr = Vector{ComplexF64}(undef, 0)
    j_idx_arr_contr = Int64[]
    tri_vals = [hcat(cos.(ang), sin.(ang)) for ang in angles]
    for (i, j) in enumerate(j_idx_arr)
        coeff = Complex(0)
        for k in 1:size(angle_history_arr[i])[1]
            tri_string = trigonometric_history_arr[i][k,:]
            phase_factor = im^(sum(tri_string))
            angle_string = angle_history_arr[i][k,:] .+ 1
            tri_string .+= 1
            coeff += prod([tri_vals[m][angle_string[m], tri_string[m]] for m in 1:M]) * phase_factor
        end
        if !isapprox(abs2(coeff), 0.0, atol=1e-16)
            push!(coeff_arr, coeff)
            push!(j_idx_arr_contr, j)
        end
    end
    return j_idx_arr_contr, coeff_arr
end


#= function symbolic_final_state_projection(M, j1, j2)
    M = convert(Int64, M)::Int64 # number of roundtrips
    j1 = convert(Int64, j1)::Int64 # first photon final state ket index
    j2 = convert(Int64, j2)::Int64 # second photon final state ket index
    @argcheck M > 0
    l1, c1 = j2lc(j1)
    l2, c2 = j2lc(j2)
    j_idx_arr_fs_1, trigonometric_history_arr_fs_1, angle_history_arr_fs_1 = symbolic_final_state_projection_sp(M, j1)
    j_idx_arr_fs_2, trigonometric_history_arr_fs_2, angle_history_arr_fs_2 = symbolic_final_state_projection_sp(M, j2)
end =#