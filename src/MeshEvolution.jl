export lcmk2j, j2lcmk, shift_timebins, beam_splitter_operator, coin_operator




function lcmk2j(N, l, c, m, k) # transforms the one index notation of basis elements to the two index notation
    l = convert(Int64, l)
    c = convert(Int64, c)
    m = convert(Int64, m)
    k = convert(Int64, k)
	n_loops = 2
    return (l * n_loops + c) * N * n_loops + m * n_loops + k + 1
end

function j2lcmk(N, j)
    j = convert(Int64,j)
	n_loops = 2
    n_loops2 = 4
    l, j = divrem(j, n_loops2 * N)
	#l = j ÷ (n_loops2 * N)
	#j -= l * n_loops2 * N
    c, j = divrem(j, n_loops * N)
	#c = j ÷ (n_loops * N)
	#j -= c * n_loops * N
    m, k = divrem(j, n_loops)
	#m = j ÷ n_loops
	#j -= m * n_loops
	#k = j
	return l, c, m, k
end

function lm2j(N, l, m): # transforms the one index notation of basis elements to the two index notation
    return l*N+m
end

function j2lm(N, j): # inverse of lm2j
    l,m = divrem(j,N)
    return l,m
end

#= function shift_timebins(state_vec::Vector)
    n_loops = 2
    new_vec = Vector{ComplexF64}(undef, length(state_vec)+n_loops)
    new_vec[2] = 0
    new_vec[end-1] = 0
    new_vec[1:2:end-3] = state_vec[1:2:end]
    new_vec[4:2:end] = state_vec[2:2:end]
    return new_vec
end =#

function shift_timebins(state_vec::Vector)
    n_loops = 2
    N = Int64(sqrt(length(state_vec)/(n_loops)))
    new_vec = zeros{ComplexF64}(length(state_vec)+n_loops)
    for j in 1:length(state_vec)
        l,c,m,k = j2lcmk(N,j)
        shifted_j = lcmk2j(N, l+c, c, m+k, k)
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

function correlated_timebin_state(coeffs::Vec)
N = length(coeffs)


function insert_initial_state()
end
