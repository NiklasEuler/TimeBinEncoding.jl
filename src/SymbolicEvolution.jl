export symbolic_generic_ket_evolution

function symbolic_generic_ket_evolution(M)
    trigonometric_history_arr = fill(Matrix{Int64}(undef, 0,M), 2*M)
    trigonometric_history_arr[1] = fill(-1, 1, M)
    angle_history_arr = copy(trigonometric_history_arr)
    for m in 1:M
        for j in 1:n_loops:m*n_loops
            l,c = j2lc(j)
            trigonometric_history_short = copy(trigonometric_history_arr[j]) # copy of trigonometric history of short loop bin to facilitate coupling
            trigonometric_history_long = copy(trigonometric_history_arr[j+1]) # copy of trigonometric history of long loop bin to facilitate coupling
            trigonometric_history_arr[j][:,m] .= 0 # cos to stay in short loop
            trigonometric_history_arr[j+1][:,m] .= 0 # cos to stay in long loop
            trigonometric_history_short[:,m] .= 1 # sin to switch to long loop
            trigonometric_history_long[:,m] .= 1 # sin to switch to short loop
            trigonometric_history_arr[j] = vcat(trigonometric_history_arr[j], trigonometric_history_long)
            trigonometric_history_arr[j+1] = vcat(trigonometric_history_arr[j+1], trigonometric_history_short)

            angle_history_arr[j][:,m] .= l# l'th time bin in short loop
            angle_history_arr[j+1][:,m] .= l # l'th time bin in short loop
            angle_history_short = copy(angle_history_arr[j]) # copy of angle history of short loop bin to facilitate coupling
            angle_history_long = copy(angle_history_arr[j+1]) # copy of angle history of long loop bin to facilitate coupling
            angle_history_arr[j] = vcat(angle_history_arr[j], angle_history_long)
            angle_history_arr[j+1] = vcat(angle_history_arr[j+1], angle_history_short)
        end
        
    end
    return trigonometric_history_arr, angle_history_arr
end