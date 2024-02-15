export lcmk2j, j2lcmk, lm2j, j2lm, lc2j, j2lc


function lcmk2j(N, l, c, m, k) # transforms the one index notation of basis elements to the two index notation
    N = convert(Int64, N)::Int64
    l = convert(Int64, l)::Int64
    c = convert(Int64, c)::Int64
    m = convert(Int64, m)::Int64
    k = convert(Int64, k)::Int64

    @argcheck N > 0
    @argcheck l ≥ 0
    @argcheck l < N
    @argcheck m ≥ 0
    @argcheck m < N
    @argcheck c ∈ [0,1]
    @argcheck k ∈ [0,1]

    return (l * n_loops + c) * N * n_loops + m * n_loops + k + 1
end

function j2lcmk(N, j)
    N = convert(Int64, N)::Int64
    j = convert(Int64, j)::Int64

    @argcheck N > 0
    @argcheck j > 0
    @argcheck j ≤ 4*N^2

    l, j = divrem(j-1, n_loops2 * N)
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

function lm2j(N, l, m) # transforms the one index notation of basis elements to the two index notation
    N = convert(Int64, N)::Int64
    l = convert(Int64, l)::Int64
    m = convert(Int64, m)::Int64

    @argcheck N > 0
    @argcheck l ≥ 0
    @argcheck l < N
    @argcheck m ≥ 0
    @argcheck m < N

    return l * N + m + 1
end

function j2lm(N, j) # inverse of lm2j
    N = convert(Int64, N)::Int64
    j = convert(Int64, j)::Int64

    @argcheck N > 0
    @argcheck j > 0
    @argcheck j ≤ N^2

    l, m = divrem(j-1, N)
    return l, m
end

function lc2j(l, c) 
    l = convert(Int64, l)::Int64
    c = convert(Int64, c)::Int64

    @argcheck l ≥ 0
    @argcheck c ∈ [0,1]
    return l * 2 + c + 1
end

function j2lc(j) # inverse of lm2j
    j = convert(Int64, j)::Int64

    @argcheck j > 0
    l,c = divrem(j-1,2)
    return l, c
end
