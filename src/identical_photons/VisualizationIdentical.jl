export visual_meas_coh_map_identical

function visual_meas_coh_map_identical end

function visual_meas_coh_map_identical(
    j_out::Int64, angles, extract_diagonal=true
)
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    j1_arr, j2_arr, weights =  explicit_fs_coherence_map_identical(j_out, angles)

    l1, c1, m1, k1, l2, c2, m2, k2 = j_super2lcmk_identical(N + M, j_out)
    c1_str = c1 == 0 ? "S" : "L"
    c2_str = c2 == 0 ? "S" : "L"
    k1_str = k1 == 0 ? "S" : "L"
    k2_str = k2 == 0 ? "S" : "L"

	println("⟨", l1, " ", c1_str, " ", m1, " ", k1_str, " ", l2, " ", c2_str, " ", m2, " ",
        k2_str, "|(SC)^M ρ (C^†S^†)^M) |", l1, " ", c1_str, " ", m1, " ", k1_str, " ", l2, " ",
        c2_str, " ", m2, " ", k2_str, "⟩ ="
    )
    _visualize_coherence_identical(N, j1_arr, j2_arr, weights, extract_diagonal)
    return nothing
end

function visual_meas_coh_map_identical(
    j_out_arr::Vector{Int64},
    angles,
    extract_diagonal=true,
    projector_weights=ones(Float64, length(j_out_arr))
)
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map_identical(j_out_arr, angles, projector_weights)

    for (j_idx, j_out) in enumerate(j_out_arr)
        l1, c1, m1, k1, l2, c2, m2, k2 = j_super2lcmk_identical(N + M, j_out)
        c1_str = c1 == 0 ? "S" : "L"
        c2_str = c2 == 0 ? "S" : "L"
        k1_str = k1 == 0 ? "S" : "L"
        k2_str = k2 == 0 ? "S" : "L"
        sign_str = sign(projector_weights[j_idx]) == 1 ? "+" : "-"
        prefactor_str = sign_str * " " * string(abs(projector_weights[j_idx]))

        println(prefactor_str, " ⟨", l1, " ", c1_str, " ", m1, " ", k1_str, " ", l2, " ",
            c2_str, " ", m2, " ", k2_str, "|(SC)^M ρ (C^†S^†)^M) |", l1, " ", c1_str, " ",
            m1, " ", k1_str, " ", l2," ", c2_str, " ", m2, " ", k2_str, "⟩ ="
        )
    end
    println("=")
    _visualize_coherence_identical(N, j1_arr, j2_arr, weights, extract_diagonal)
    return nothing
end

function _visualize_coherence_identical(N, j1_arr, j2_arr, weights, extract_diagonal=true)
    contr_j_idxs = correlated_short_bins_idxs_identical(N)
    extractable_correlated_coherences = []
    display_weights = round.(Real.(weights), digits=5)

	for i in eachindex(j1_arr)
        j1 = j1_arr[i]
        j2 = j2_arr[i]

		l1_bra, c1_bra, m1_bra, k1_bra, l2_bra, c2_bra, m2_bra, k2_bra =
            j_super2lcmk_identical(N, j1)
		l1_ket, c1_ket, m1_ket, k1_ket, l2_ket, c2_ket, m2_ket, k2_ket =
            j_super2lcmk_identical(N, j2)

        if(j1 in contr_j_idxs && j2 in contr_j_idxs && (extract_diagonal || j1 ≠ j2))
            push!(extractable_correlated_coherences, i)
        end
		println("+ ρ_[", l1_bra, " ", m1_bra, " ", l2_bra, " ", m2_bra, "]^[", l1_ket, " ",
            m1_ket, " ", l2_ket, " ", m2_ket, "] ⋅ ", display_weights[i]
        )
	end
    println("Useful extractable coherences:")
    for i in extractable_correlated_coherences
        j1 = j1_arr[i]
        j2 = j2_arr[i]
		l1_bra, c1_bra, m1_bra, k1_bra, l2_bra, c2_bra, m2_bra, k2_bra =
            j_super2lcmk_identical(N, j1)
		l1_ket, c1_ket, m1_ket, k1_ket, l2_ket, c2_ket, m2_ket, k2_ket =
            j_super2lcmk_identical(N, j2)
        println("+ ρ_[", l1_bra, " ", m1_bra, " ", l2_bra, " ", m2_bra, "]^[", l1_ket, " ",
            m1_ket, " ", l2_ket, " ", m2_ket, "] ⋅ ", display_weights[i]
        )
   end

    return nothing
end
