#= export symbolic_ket_evolution

function symbolic_generic_ket_evolution(M)
    state_indicator = [1.0,0.0]
    coupling_history_arr = zeros(M, 2*M)
    timebin_history_arr = zero(coupling_history_arr)
    for m in 1:M
        for j in 1:2*M
            l,c = j2lc(j)

    end

end =#