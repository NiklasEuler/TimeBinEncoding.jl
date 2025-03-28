export lcmk2j, j2lcmk, lm2j, j2lm, lc2j, j2lc
#export lcmk2j_identical, j2lcmk_identical
#export lcmk2j_super_identical, j_super2lcmk_identical, correlated_short_bins_idxs_identical
#export correlated_short_bins_tuples_identical
export correlated_short_bins_idxs, correlated_short_bins_tuples, indices2tuples

"""
    lcmk2j(N, l, c, m, k)

Convert the parameters of the two-photon state `|lcmk⟩` to its corresponding
basis index `j`.

# Arguments
- `N`: number of time bins in the basis
- `l`: time bin index of the signal photon, `l` ∈ {0, 1,…}, where `l` indicates the total
    number of long roundtrips taken.
- `c`: loop index of the signal photon, `c` ∈ {0, 1}, where `c==0` means the short loop and
    `c== 1` means the long loop.
- `m`: time bin index of the idler photon, `m` ∈ {0, 1,…}, with the same encoding as above.
- `k`: loop index of the idler photon, `k` ∈ {0, 1}, with the same encoding as above.

See also [`j2lcmk`](@ref), [`lm2j`](@ref), [`j2lm`](@ref), [`lc2j`](@ref), [`j2lc`](@ref).
"""
function lcmk2j(N, l, c, m, k)
    # transforms the one index notation of basis elements to the two index notation

    N, l, c, m, k = _lcmk2j_input_sanity_check(N, l, c, m, k)

    return (l * N_LOOPS + c) * N * N_LOOPS + m * N_LOOPS + k + 1
end

function _lcmk2j_input_sanity_check(N, l, c, m, k)
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
    @argcheck c in [0, 1]
    @argcheck k in [0, 1]

    return N, l, c, m, k
end


"""
    j2lcmk(N, j)

Inverse function of `lcmk2j`.

See also [`lcmk2j`](@ref), [`lm2j`](@ref), [`j2lm`](@ref), [`lc2j`](@ref), [`j2lc`](@ref).
"""
function j2lcmk(N, j)

    N, j = _j2lcmk_input_sanity_check(N, j)

    l, j = divrem(j - 1, N_LOOPS2 * N)
    c, j = divrem(j, N_LOOPS * N)
    m, k = divrem(j, N_LOOPS)
	return l, c, m, k
end

function _j2lcmk_input_sanity_check(N, j)
    N = convert(Int64, N)::Int64
    j = convert(Int64, j)::Int64

    @argcheck N > 0
    @argcheck j > 0
    @argcheck j ≤ N_LOOPS2 * N ^ 2

    return N, j
end

"""
    lm2j(N, l, m)

Convert two indices `l`,`m` ∈ {0,…, N - 1} to a joint index `j`.

See also [`j2lm`](@ref), [lcmk2j`](@ref), [`j2lcmk`](@ref), [`lc2j`](@ref), [`j2lc`](@ref).
"""
function lm2j(N, l, m)
    # transforms the one index notation of basis elements to the two index notation

    N, l, m = _lm2j_input_sanity_check(N, l, m)

    return l * N + m + 1
end

function _lm2j_input_sanity_check(N, l, m)
    N = convert(Int64, N)::Int64
    l = convert(Int64, l)::Int64
    m = convert(Int64, m)::Int64

    @argcheck N > 0
    @argcheck l ≥ 0
    @argcheck l < N
    @argcheck m ≥ 0
    @argcheck m < N

    return N, l, m
end

"""
    j2lm(N, j)

Inverse function of `lm2j`.

See also [`j2lm`](@ref), [`lcmk2j`](@ref), [`j2lcmk`](@ref), [`lc2j`](@ref), [`j2lc`](@ref).
"""
function j2lm(N, j) # inverse of lm2j

    N, j = _j2lm_input_sanity_check(N, j)

    l, m = divrem(j - 1, N)
    return l, m
end

function _j2lm_input_sanity_check(N, j)
    N = convert(Int64, N)::Int64
    j = convert(Int64, j)::Int64

    @argcheck N > 0
    @argcheck j > 0
    @argcheck j ≤ N^2

    return N, j
end

"""
    lc2j(l, c)

Convert the parameters of the single-photon state `|lc⟩` to its corresponding
basis index `j`.

# Arguments
- `l`: time bin index of the photon, `l ∈ {0, 1, …}`, where `l` indicates the total number
    of long roundtrips taken.
- `c`: loop index of the photon, `c ∈ {0, 1}`, where `c == 0` means the short loop and
    `c == 1` means the long loop.

See also [`j2lc`](@ref), [`lcmk2j`](@ref), [`j2lcmk`](@ref), [`lm2j`](@ref), [`j2lm`](@ref).

"""
function lc2j(l, c)

    l,c = _lc2j_input_sanity_check(l, c)

    return l * 2 + c + 1
end

function _lc2j_input_sanity_check(l, c)
    l = convert(Int64, l)::Int64
    c = convert(Int64, c)::Int64

    @argcheck l ≥ 0
    @argcheck c in [0, 1]

    return l, c
end

"""
    j2lc(j)

Inverse function of `lc2j`.

See also [`lc2j`](@ref), [`lcmk2j`](@ref), [`j2lcmk`](@ref), [`lm2j`](@ref), [`j2lm`](@ref).
"""
function j2lc(j) # inverse of lm2j
    j = convert(Int64, j)::Int64

    @argcheck j > 0

    l, c  = divrem(j - 1, 2)
    return l, c
end

"""
    correlated_short_bins_idxs(N)

Compute the joint indices `j` of the `|lcmk⟩` basis corresponding to `|i0i0⟩` for
i ∈ {0,…, N - 1}.

See also [`correlated_short_bins_tuples`](@ref), [`indices2tuples`](@ref).
"""
function correlated_short_bins_idxs(N)
    N = convert(Int64, N)::Int64
    contr_j_idxs = [lcmk2j(N, i, 0, i, 0) for i in 0:N - 1]
    return contr_j_idxs
end



"""
    correlated_short_bins_tuples(N; extract_diagonal=false)

Convert the indices of correlated short bins to tuples.

# Arguments
- `N`: The number of time bins.
- `extract_diagonal`: (Optional) A boolean indicating whether to include the diagonal
    elements. Default is `false`.

# Returns
- `Vector{Tuple{Int64, Int64}}`: A list of tuples representing pairs of indices.

See also [`correlated_short_bins_idxs`](@ref), [`indices2tuples`](@ref).
"""
function correlated_short_bins_tuples(N; extract_diagonal=false)
    contr_j_idxs = correlated_short_bins_idxs(N)
    return indices2tuples(contr_j_idxs; extract_diagonal=extract_diagonal)
end

"""
    indices2tuples(indices; extract_diagonal=false)

Converts a list of indices into a list of tuples. Each tuple represents a pair of indices
(i, j) where i and j are elements from the input list `indices`. By default, the function
includes all possible pairs of indices but excludes the diagonal elements (i, i). They can
also be included by setting `extract_diagonal` to `true`.

# Arguments
- `indices`: A list of indices.
- `extract_diagonal`: (Optional) A boolean value indicating whether to include the diagonal
    elements in the output. Default is `false`.

# Returns
- `Vector{Tuple{Int64, Int64}}`: A list of tuples representing pairs of indices.

See also [`correlated_short_bins_tuples`](@ref), [`correlated_short_bins_idxs`](@ref).
"""
function indices2tuples(indices; extract_diagonal=false)
    return [(i,j) for i in indices for j in indices if i ≠ j || extract_diagonal]
end
