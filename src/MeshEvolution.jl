export lcmk2j, lcmk2j, shift_timebins, beam_splitter_operator, coin_operator




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
	l = j ÷ (n_loops2 * N)
	j -= l * n_loops2 * N
	c = j ÷ (n_loops * N)
	j -= c * n_loops * N
	m = j ÷ n_loops
	j -= m * n_loops
	k = j
	return l, c, m, k
end

function shift_timebins(state_vec::Vector)
    new_vec = Vector{ComplexF64}(undef, length(state_vec)+4)
    new_vec[2] = 0
    new_vec[end-1] = 0
    new_vec[1:2:end-3] = state_vec[1:2:end]
    new_vec[2:2:end-1] = state_vec[2:2:end]
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
    return blockdiag(matrices...)
end

#function insert_initial_state()