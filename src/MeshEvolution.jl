export shift_timebins, beam_splitter_operator, coin_operator, mesh_evolution
export correlated_timebin_state, insert_initial_state
export shift_timebins_single_photon

function shift_timebins_single_photon(state_vec::Vector)
    new_vec = Vector{ComplexF64}(undef, length(state_vec)+n_loops)
    new_vec[2] = 0
    new_vec[end-1] = 0
    new_vec[1:2:end-3] = state_vec[1:2:end]
    new_vec[4:2:end] = state_vec[2:2:end]
    return new_vec
end

function shift_timebins(state_vec::Vector)
    N = Int64(sqrt(length(state_vec)/(n_loops2)))
    new_vec = zeros(ComplexF64, ((N+1)*n_loops)^2)
    for j in 1:length(state_vec)
        l,c,m,k = j2lcmk(N,j)
        shifted_j = lcmk2j(N+1, l+c, c, m+k, k) # adapted system has one more time bin, so we need to put N+1
        new_vec[shifted_j] = state_vec[j]
    end
    return new_vec
end

function beam_splitter_operator(θ)
    θ = convert(Float64,θ)
    cs = cos(θ)
    sn = im*sin(θ)
    cols = [1,1,2,2]
    rows = [1,2,1,2]
    vals = [cs, sn, sn, cs]
   return sparse(cols,rows,vals)
end

function coin_operator(angles::Vector)
    real_angles = convert(Vector{Float64}, angles)
    matrices = [beam_splitter_operator(θ) for θ in real_angles]
    single_photon_coin_operator = blockdiag(matrices...)
    tensor_coin_operator = kron(single_photon_coin_operator,single_photon_coin_operator)
    return tensor_coin_operator
end

function correlated_timebin_state(wf_coeffs::Vector)#todo: normalize
    wf_coeffs = convert(Vector{ComplexF64}, wf_coeffs)
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
    state = copy(ψ_init)
    for i in eachindex(angles)
        coin_op = coin_operator(angles[i])
        state = coin_op * state
        state = shift_timebins(state)
    end
    return state
end