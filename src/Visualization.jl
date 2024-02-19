export visualize_symbolic_ket_evolution_sp, visualize_symbolic_final_state_projection_sp, visualize_measurement_coherence_map


"""
    visualize_symbolic_ket_evolution_sp(M, l_init)

TBW
"""
function visualize_symbolic_ket_evolution_sp(M, l_init)
    j_idx_arr, trigonometric_history_arr, angle_history_arr = symbolic_ket_evolution_sp(M, l_init)
    println("The initial state |$l_init,S⟩ after $M roundtrips has evolved into the following contributions:")
    println("")
    println("(SC)^{$M} |$l_init,S⟩ =")
    cs = ["S", "L"]
    for (i, j) in enumerate(j_idx_arr)
        l,c = j2lc(j)
        print("+ |$l,$(cs[c+1])⟩⋅[")
        trigonometric_string_formatter(trigonometric_history_arr[i], angle_history_arr[i])
        println(" ]")
    end
end

"""
    visualize_symbolic_final_state_projection_sp(M, l_fs, c_fs)

TBW
"""
function visualize_symbolic_final_state_projection_sp(M, l_fs, c_fs)
    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs = symbolic_final_state_projection_sp(M, l_fs, c_fs)
    #println("The initial state |$l_init,S⟩ after $M roundtrips has evolved into the following contributions:")
    #println("")
    cs = ["S", "L"]
    println("|$l_fs,$(cs[c_fs+1])⟩⟨$l_fs,$(cs[c_fs+1])| (SC)^{$M} |Ψ_init⟩ = |$l_fs,$(cs[c_fs+1])⟩⋅{")
    for (i, j) in enumerate(j_idx_arr_fs)
        l_init,c_init = j2lc(j)
        print("+ c_($l_init,$(cs[c_init+1]))⋅[")
        trigonometric_string_formatter(trigonometric_history_arr_fs[i], angle_history_arr_fs[i])
        println(" ]")
    end
    println("}")
end

"""
    trigonometric_string_formatter(trigonometric_history, angle_history)

TBW
"""
function trigonometric_string_formatter(trigonometric_history, angle_history)
    n_tri_strings = size(trigonometric_history)[1]
    tri_replacement_str = ["cos(θ","sin(θ"]
    phase_replacement_str = [" + "," + i ⋅ "," - "," - i ⋅ "]
    for k in 1:n_tri_strings
        tri_string_bin = trigonometric_history[k,:]
        angle_string_bin = angle_history[k,:]
        tri_string_formatted = join([join([tri_replacement_str[1+tri_string_bin[i]],"_$(angle_string_bin[i])^$i)"]) for i in eachindex(tri_string_bin)])
        phase_tri_string = join([phase_replacement_str[mod(sum(tri_string_bin),4)+1],tri_string_formatted])
        print(phase_tri_string)
    end
end

function visualize_measurement_coherence_map(j_out, angles)
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    j1_arr, j2_arr, weights =  explicit_measurement_coherence_map(j_out, angles)
    display_weights = round.(Real.(weights), digits=5)
    l_out,c_out,m_out,k_out = j2lcmk(N+M,j_out)
	println("⟨",l_out," ",c_out," ",m_out," ",k_out,"|(SC)^M ρ (C^†S^†)^M) |",l_out," ",c_out," ",m_out," ",k_out,"⟩ =")
	for i in eachindex(j1_arr)
		l1, c1, m1, k1 = j2lcmk(N,j1_arr[i])
		l2, c2, m2, k2 = j2lcmk(N,j2_arr[i])
		println("+ ρ_[",l1," ",m1,"]^[",l2," ",m2,"] ⋅ ",display_weights[i])
	end
end