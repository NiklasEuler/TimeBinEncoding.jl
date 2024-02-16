export visualize_symbolic_ket_evolution_sp, visualize_symbolic_final_state_projection_sp


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

