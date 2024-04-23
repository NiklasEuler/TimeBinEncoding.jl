export visualize_symbolic_ket_evolution_sp, visualize_symbolic_fs_projection_sp,
    visualize_measurement_coherence_map, visualize_combined_measurement_coherence_map
export visualize_measurement_coherence_map_identical

"""
    visualize_symbolic_ket_evolution_sp(M, l_init)

TBW
"""
function visualize_symbolic_ket_evolution_sp(M, l_init)
    j_idx_arr, trigonometric_history_arr, angle_history_arr =
        symbolic_ket_evolution_sp(M, l_init)
    println("The initial state |$l_init,S⟩ after $(M) "*
        "roundtrips has evolved into the following contributions:")
    println("")
    println("(SC)^{$M} |$l_init,S⟩ =")
    cs = ["S", "L"]
    for (i, j) in enumerate(j_idx_arr)
        l, c  = j2lc(j)
        print("+ |$l,$(cs[c + 1])⟩⋅[")
        trigonometric_string_formatter(trigonometric_history_arr[i], angle_history_arr[i])
        println("]")
   end

    return nothing
end

"""
    visualize_symbolic_fs_projection_sp(M, l_fs, c_fs)

TBW
"""
function visualize_symbolic_fs_projection_sp(M, l_fs, c_fs)
    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs =
        symbolic_fs_projection_sp(M, l_fs, c_fs)
    cs = ["S", "L"]
    println("|$l_fs,$(cs[c_fs + 1])⟩⟨$l_fs,$(cs[c_fs + 1])| (SC)^{$M} |Ψ_init⟩ = |$(l_fs)"*
        ",$(cs[c_fs + 1])⟩⋅{")
    for (i, j) in enumerate(j_idx_arr_fs)
        l_init, c_init = j2lc(j)
        print("+ c_($l_init,$(cs[c_init + 1]))⋅[")
        trigonometric_string_formatter(trigonometric_history_arr_fs[i],
            angle_history_arr_fs[i])
        println("]")
    end
    println("}")
    return nothing
end

"""
    trigonometric_string_formatter(trigonometric_history, angle_history)

TBW
"""
function trigonometric_string_formatter(trigonometric_history, angle_history)
    n_tri_strings = size(trigonometric_history)[1]
    tri_replacement_str = ["cos(θ", "sin(θ"]
    phase_replacement_str = [" + ", " + i ⋅ ", " - ", " - i ⋅ "]
    for k in 1:n_tri_strings
        tri_string_bin = trigonometric_history[k, :]
        angle_string_bin = angle_history[k, :]
        tri_string_formatted = join([join([tri_replacement_str[1 + tri_string_bin[i]],
            "_$(angle_string_bin[i])^$i)"]) for i in eachindex(tri_string_bin)])
        phase_tri_string =
            join([phase_replacement_str[mod(sum(tri_string_bin), 4) + 1],
                tri_string_formatted]
            )
        print(phase_tri_string)
   end

    return nothing
end

function visualize_measurement_coherence_map end

function visualize_measurement_coherence_map(
    j_out::Int64, angles, extract_diagonal=true, off_l=0, off_m=0)
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    j1_arr, j2_arr, weights =  explicit_fs_coherence_map(j_out, angles)

    if(off_l != 0 || off_m != 0)
        println("with offsets: off_l = ", off_l, ", off_m = ", off_m)
    end

    l_out, c_out, m_out, k_out = j2lcmk(N + M, j_out)
	println("⟨", l_out, " ", c_out, " ", m_out, " ", k_out, "|(SC)^M ρ (C^†S^†)^M) |",
        l_out, " ", c_out, " ", m_out, " ", k_out, "⟩ ="
    )
    _visualize_coherence(N, j1_arr, j2_arr, weights, extract_diagonal, off_l, off_m)
    return nothing
end

function visualize_measurement_coherence_map(
    j_out_arr::Vector{Int64},
    angles,
    extract_diagonal=true,
    projector_weights=ones(Float64, length(j_out_arr)),
    off_l=0,
    off_m=0
)
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map(j_out_arr, angles, projector_weights)

    if(off_l != 0 || off_m != 0)
        println("with offsets: off_l = ", off_l, ", off_m = ", off_m)
    end

    for (j_idx, j_out) in enumerate(j_out_arr)
        l_out, c_out, m_out, k_out = j2lcmk(N + M, j_out)
        c_str = c_out == 0 ? "S" : "L"
        k_str = k_out == 0 ? "S" : "L"
        sign_str = sign(projector_weights[j_idx]) == 1 ? "+" : "-"
        prefactor_str = sign_str * " " * string(abs(projector_weights[j_idx]))

        println(prefactor_str, " ⟨", l_out, " ", c_str, " ", m_out, " ", k_str,
            "|(SC)^M ρ (C^†S^†)^M) |", l_out, " ", c_str, " ", m_out, " ", k_str, "⟩"
        )
    end

    println("=")
    _visualize_coherence(N, j1_arr, j2_arr, weights, extract_diagonal, off_l, off_m)

    return nothing
end



function _visualize_coherence(
    N,
    j1_arr,
    j2_arr,
    weights,
    extract_diagonal=true,
    off_l=0,
    off_m=0
)
    contr_j_idxs = correlated_short_bins_idxs(N)
    extractable_correlated_coherences = []
    display_weights = round.(Real.(weights), digits=5)

	for i in eachindex(j1_arr)
        j1 = j1_arr[i]
        j2 = j2_arr[i]
		l1, c1, m1, k1 = j2lcmk(N, j1)
		l2, c2, m2, k2 = j2lcmk(N, j2)

        l1 -= off_l
        l2 -= off_l
        m1 -= off_m
        m2 -= off_m
        if(l1 < 0 || l2 < 0 || m1 < 0 || m2 < 0)
            continue
        end

        j1 = lcmk2j(N, l1, c1, m1, k1)
        j2 = lcmk2j(N, l2, c2, m2, k2)

        if(j1 in contr_j_idxs && j2 in contr_j_idxs && (extract_diagonal || j1 ≠ j2))
            push!(extractable_correlated_coherences, i)
        end
		println("+ ρ_[", l1, " ", m1, "]^[", l2, " ", m2, "] ⋅ ", display_weights[i])
	end
    println("Useful extractable coherences:")
    for i in extractable_correlated_coherences
        j1 = j1_arr[i]
        j2 = j2_arr[i]
		l1, c1, m1, k1 = j2lcmk(N, j1)
		l2, c2, m2, k2 = j2lcmk(N, j2)
        println("+ ρ_[", l1 - off_l, " ", m1 - off_m, "]^[", l2 - off_l, " ", m2 - off_m,
            "] ⋅ ", display_weights[i]
        )
   end

    return nothing
end

function visualize_measurement_coherence_map_identical end

function visualize_measurement_coherence_map_identical(
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

function visualize_measurement_coherence_map_identical(
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
