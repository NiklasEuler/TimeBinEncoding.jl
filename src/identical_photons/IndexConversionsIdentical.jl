export lcmk2j_identical, j2lcmk_identical
export lcmk2j_super_identical, j_super2lcmk_identical, correlated_short_bins_idxs_identical
export correlated_short_bins_tuples_identical

function correlated_short_bins_idxs_identical(N)
    contr_j_idxs = Int64[]
    N = convert(Int64, N)::Int64

    for i in 0:N - 1
        for j in i:N - 1
            push!(contr_j_idxs, lcmk2j_super_identical(N, i, 0, j, 0, i, 0, j, 0))
        end
    end

    return contr_j_idxs
end

function correlated_short_bins_tuples_identical(N; extract_diagonal=false)
    contr_j_idxs = correlated_short_bins_idxs_identical(N)
    return [(i,j) for i in contr_j_idxs for j in contr_j_idxs if i ≠ j || extract_diagonal]
end

function lcmk2j_identical(N, l, c, m, k)

    N, l, c, m, k = _lcmk2j_input_sanity_check(N, l, c, m, k)

    @argcheck m ≥ l

    if l == m && c == 1 && k == 0
        throw(ArgumentError("For identical time bins, the loop indices must be ordered due
        to indistinguishability."))
    end

    dynamic_base_width = append!([0], cumsum(N_LOOPS * N:-1:1))
    # cumulative sum of the number of basis elements for each time bin
    # reduces in size increment since m ≥ l and index skip for LS loops when m == l

    j = dynamic_base_width[l * N_LOOPS + c + 1] + (m - l) * N_LOOPS + k - c + 1
    # -c to skip one j index in case of identical time bins
    # ToDo: Need to extend documentation on this function.
    return j

end

function j2lcmk_identical(N, j)

    N, j = _j2lcmk_input_sanity_check(N, j)

    dynamic_base_width = append!([0], cumsum(N_LOOPS * N:-1:1))
    # cumulative sum of the number of basis elements for each time bin
    # reduces in size increment since m ≥ l and index skip for LS loops when m == l

    j = j - 1 # correct for 1-based indexing
    i = 1
    while j >= dynamic_base_width[i + 1]
        i += 1
    end

    l = i - 1
    j -= dynamic_base_width[i]

    l, c = divrem(l, N_LOOPS)

    j += c + l * N_LOOPS

    m, k = divrem(j, N_LOOPS)


   #=  l_temp, c_temp, m_temp, k_temp = j2lcmk(N, j)
    idx_correction = l_temp
    # for each timebin there is one less possible index due to indistinguishability

    l_temp, c_temp, m_temp, k_temp = j2lcmk(N, j + idx_correction)
    # correct for clear shift at this point

    idx_correction = l_temp # update in case of correction pushing over l_temp

    if m_temp ≥ l_temp && c_temp == 1
    # possibly one additional index shift in current m cycle
        idx_correction += 1
    end

    l, c, m, k = j2lcmk(N, j + idx_correction) # final indices =#
	return l, c, m, k
end



"""
    j_super2lcmk_identical(N, j_super)

Converts a 4-photon basis state index ´j_super´ in the |l1, c1 , m1 , k1, l2, c2 , m2 , k2>
basis to the corresponding |l1, c1 , m1 , k1> and |l2, c2 , m2 , k2> indices for the signal
and idler photon pairs, respectively.

"""
function j_super2lcmk_identical(N, j_super)
    N = convert(Int64, N)::Int64
    j_super = convert(Int64, j_super)::Int64
    # 4-photon bin index in the |l1, c1 , m1 , k1, l2, c2 , m2 , k2 basis

    d_hilbert_space = Int(N_LOOPS * N * (N_LOOPS * N + 1) / 2)
    # local hilbert space dimension for two photons
    # indistinguishable and N_LOOPS * N states per photon -> geometric series
    j1, j2 = j2lm(d_hilbert_space, j_super)
    # translate j_super to two 2-photon j indices (indistinguishable/identical )
    l1, c1 , m1 , k1  = j2lcmk_identical(N, j1 + 1) # +1 for 1-based indexing
    l2, c2 , m2 , k2  = j2lcmk_identical(N, j2 + 1) # +1 for 1-based indexing

    return l1, c1, m1, k1, l2, c2, m2, k2

end

"""
    lcmk2j_super_identical(N, l1, c1, m1, k1, l2, c2, m2, k2)

    Converts the |l1, c1 , m1 , k1> and |l2, c2 , m2 , k2> indices for the signal and idler
    photon pairs, respectively, to the corresponding 4-photon basis state index ´j_super´ in
    the |l1, c1 , m1 , k1, l2, c2 , m2 , k2> basis.

"""
function lcmk2j_super_identical(N, l1, c1, m1, k1, l2, c2, m2, k2)

    j1 = lcmk2j_identical(N, l1, c1, m1, k1)
    j2 = lcmk2j_identical(N, l2, c2, m2, k2)

    d_hilbert_space = Int(N_LOOPS * N * (N_LOOPS * N + 1) / 2)
        # local hilbert space dimension for two photons
    # indistinguishable and N_LOOPS * N states per photon -> geometric series
    j_super = lm2j(d_hilbert_space, j1 - 1, j2 - 1) # -1 for 0-based indexing

    return j_super
end
