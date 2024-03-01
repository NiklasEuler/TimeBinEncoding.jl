export lcmk2j, j2lcmk, lm2j, j2lm, lc2j, j2lc, correlated_short_bins_idxs


"""
    lcmk2j(N, l, c, m, k)

Convert the parameters of the two-photon state '|lcmk⟩' to its corresponding basis index `j`.

# Arguments
- `N`: number of time bins in the basis
- `l`: time bin index of the signal photon, `l` ∈ {0,1,…}, where `l` indicates the total number of long roundtrips taken.
- `c`: loop index of the signal photon, `c` ∈ {0,1}, where `c==0` means the short loop and `c==1` the long loop.
- `m`: time bin index of the idler photon, `m` ∈ {0,1,…}, with the same encoding as above.
- `k`: loop index of the idler photon, `k` ∈ {0,1}, with the same encoding as above.

See also `j2lcmk`, `lm2j`, `j2lm`, `lc2j`, `j2lc`.
"""
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

function correlated_short_bins_idxs(N)
    N = convert(Int64, N)::Int64
    contr_j_idxs = [lcmk2j(N,i,0,i,0) for i in 0:N-1]
    return contr_j_idxs
end