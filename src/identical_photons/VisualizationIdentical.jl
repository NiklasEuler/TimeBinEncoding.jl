export visual_meas_coh_map_identical
export visual_meas_coh_map_combined_identical

"""
    visual_meas_coh_map_identical(
        j_out, angles, phases=ones(length(angles[1])); extract_diagonal=true
    )

Visualize the contributions of the initial state coherences to a final-state projector given
by `j_out` in the `|l1, c1, m1, k1, l2, c2, m2, k2⟩` basis after evolution through beam-
splitter `angles`. Optionally, initial-state phases can be specified by `phases`.
At the end, a list of all contributing coherences to the fidelity of a typical reference
state is printed. Here, by default, the diagonal elements are also included, but this can be
changed by setting `extract_diagonal=false`.

See also [`visual_meas_coh_map`](@ref), [`visual_meas_coh_map_combined_identical`](@ref).

"""
function visual_meas_coh_map_identical end

function visual_meas_coh_map_identical(
    j_out::Int64, angles, phases=ones(length(angles[1])); extract_diagonal=true
)
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    projector_weight = 1 # default projector weight
    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map_identical(j_out, angles, projector_weight, phases)

    l1, c1, m1, k1, l2, c2, m2, k2 = j_super2lcmk_identical(N + M, j_out)
    c1_str = c1 == 0 ? "S" : "L"
    c2_str = c2 == 0 ? "S" : "L"
    k1_str = k1 == 0 ? "S" : "L"
    k2_str = k2 == 0 ? "S" : "L"

	println("⟨", l1, " ", c1_str, " ", m1, " ", k1_str, " ", l2, " ", c2_str, " ", m2, " ",
        k2_str, "|(SC)^M ρ (C^†S^†)^M) |", l1, " ", c1_str, " ", m1, " ", k1_str, " ", l2,
        " ", c2_str, " ", m2, " ", k2_str, "⟩ ="
    )
    _visual_coh_identical(N, j1_arr, j2_arr, weights, extract_diagonal, phases)
    return nothing
end

function visual_meas_coh_map_identical(
    j_out_arr::Vector{Int64},
    angles,
    projector_weights=ones(Float64, length(j_out_arr)),
    phases=ones(Float64, length(angles[1]))::Vector;
    extract_diagonal=true,
)
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map_identical(j_out_arr, angles, projector_weights, phases)

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
    _visual_coh_identical(N, j1_arr, j2_arr, weights, extract_diagonal, phases)
    return nothing
end

function _visual_coh_identical(
    N, j1_arr, j2_arr, weights, extract_diagonal=true, phases=ones(Float64, N)
)
    contr_j_idxs = correlated_short_bins_idxs_identical(N)
    extractable_correlated_coherences = []

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
        val = round(weights[i], digits=5)
        val = try
            convert(Float64, val)
        catch
            val
        end
		println("+ ρ_[", l1_bra, " ", m1_bra, " ", l2_bra, " ", m2_bra, "]^[", l1_ket, " ",
            m1_ket, " ", l2_ket, " ", m2_ket, "] ⋅ ", val
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
        val = round(weights[i], digits=5)
        val = try
            convert(Float64, val)
        catch
            val
        end
        println("+ ρ_[", l1_bra, " ", m1_bra, " ", l2_bra, " ", m2_bra, "]^[", l1_ket, " ",
            m1_ket, " ", l2_ket, " ", m2_ket, "] ⋅ ", val
        )
   end

    return nothing
end


"""
    visual_meas_coh_map_combined_identical(
        N, combined_weights, contr_j_idxs; extract_diagonal=false
    )

Visualize the initial-state contributions to a combination of different measurements. The
`combined_weights` are the weights of the combined measurements of `N` initial time bins.
The `contr_j_idxs` are the indices of the contributing coherences, which are printed at the
end. By default, the diagonal elements are not included, but this can be changed by setting
`extract_diagonal=true`.

See also [`visual_meas_coh_map_identical`](@ref), [`visual_meas_coh_map`](@ref).

"""
function visual_meas_coh_map_combined_identical(
    N, combined_weights, contr_j_idxs; extract_diagonal=false
)
	d_full_hs_bl = (N * (2 * N + 1)) ^ 2
	println("All remaining coherences:")
	for j_comb in eachindex(combined_weights)
		if !(isapprox(combined_weights[j_comb], 0, atol=1e-10))
            # only print non-zero coherences
			j1, j2 = j2lm(d_full_hs_bl, j_comb)
			j1 += 1
			j2 += 1
	 		l1_bra, c1_bra, m1_bra, k1_bra, l2_bra, c2_bra, m2_bra, k2_bra =
	            j_super2lcmk_identical(N, j1)
			l1_ket, c1_ket, m1_ket, k1_ket, l2_ket, c2_ket, m2_ket, k2_ket =
	            j_super2lcmk_identical(N, j2)
            val = round(combined_weights[j_comb], digits=5)
            val = try
                convert(Float64, val)
            catch
                val
            end
	        println("+ ρ_[", l1_bra, " ", m1_bra, " ", l2_bra, " ", m2_bra, "]^[", l1_ket,
                " ", m1_ket, " ", l2_ket, " ", m2_ket, "] ⋅ ", val
            )
	   end
	end
	println("\n All contributing coherences:")
	for j_contr1 in contr_j_idxs, j_contr2 in contr_j_idxs
		j_contr = lm2j(d_full_hs_bl,  j_contr1 - 1, j_contr2 - 1)
		if !(isapprox(combined_weights[j_contr], 0, atol=1e-10)) && (j_contr1 != j_contr2 ||
            extract_diagonal
        )
			l1_bra, c1_bra, m1_bra, k1_bra, l2_bra, c2_bra, m2_bra, k2_bra =
	            j_super2lcmk_identical(N, j_contr1)
			l1_ket, c1_ket, m1_ket, k1_ket, l2_ket, c2_ket, m2_ket, k2_ket =
	            j_super2lcmk_identical(N, j_contr2)
            val = round(combined_weights[j_contr], digits=5)
            val = try
                convert(Float64, val)
            catch
                val
            end
            println("+ ρ_[", l1_bra, " ", m1_bra, " ", l2_bra, " ", m2_bra, "]^[", l1_ket,
                " ", m1_ket, " ", l2_ket, " ", m2_ket, "] ⋅ ", val
            )
	        #println("+ ρ_[", l1_bra, " ", m1_bra, " ", l2_bra, " ", m2_bra, "]^[", l1_ket,
            #    " ", m1_ket, " ", l2_ket, " ", m2_ket, "] ⋅ ", val
            #)
	    end
	end
end
