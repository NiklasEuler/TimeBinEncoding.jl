export symbolic_generic_ket_evolution_sp, symbolic_ket_evolution_sp, visualize_symbolic_ket_evolution_sp

function symbolic_generic_ket_evolution_sp(M)
    M = convert(Int64, M)::Int64 # number of roundtrips
    @argcheck M > 0

    trigonometric_history_arr = fill(Matrix{Int64}(undef, 0, M), n_loops * (M + 1))
    trigonometric_history_arr[1] = fill(-1, 1, M)
    angle_history_arr = copy(trigonometric_history_arr)
    j_idx_arr = [1,collect(1+n_loops:n_loops * M)..., n_loops * (M + 1)]
    for m in 1:M
        symbolic_generic_ket_coin_sp(m, trigonometric_history_arr, angle_history_arr)
        symbolic_generic_ket_shift_sp(M, m, trigonometric_history_arr, angle_history_arr)
    end
    deleteat!(trigonometric_history_arr,[n_loops,length(trigonometric_history_arr)-1]) #  remove undef matrices for the two empty bins
    deleteat!(angle_history_arr,[n_loops,length(angle_history_arr)-1]) #  remove undef matrices for the two empty bins

    return j_idx_arr, trigonometric_history_arr, angle_history_arr
end

function symbolic_generic_ket_coin_sp(m, trigonometric_history_arr, angle_history_arr)
    m = convert(Int64, m)::Int64 # current round trip index
    @argcheck m > 0

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

function symbolic_generic_ket_shift_sp(M, m, trigonometric_history_arr, angle_history_arr)
    M = convert(Int64, M)::Int64 # number of roundtrips
    m = convert(Int64, m)::Int64 # current round trip index
    @argcheck m > 0
    @argcheck m ≤ M
    @argcheck M > 0

    for j in n_loops*(M+1):-n_loops:2*n_loops
        trigonometric_history_arr[j] = trigonometric_history_arr[j-n_loops]
        angle_history_arr[j] = angle_history_arr[j-n_loops]
    end
    trigonometric_history_arr[n_loops] = Matrix{Int64}(undef, 0, M) # reset 0L bin after shift
    angle_history_arr[n_loops] = Matrix{Int64}(undef, 0, M) # reset 0L bin after shift
end

function symbolic_ket_evolution_sp(M, l)
    M = convert(Int64, M)::Int64 # number of roundtrips
    l = convert(Int64, l)::Int64 # initial time bin index
    @argcheck M > 0
    @argcheck l ≥ 0

    j_idx_arr, trigonometric_history_arr, angle_history_arr = symbolic_generic_ket_evolution_sp(M)
    if l ≠ 0
        for i in eachindex(angle_history_arr)
            angle_history_arr[i] .+= l
        end
        j_idx_arr .+= n_loops*l
    end
    return j_idx_arr, trigonometric_history_arr, angle_history_arr
end

function visualize_symbolic_ket_evolution_sp(M, l)
    j_idx_arr, trigonometric_history_arr, angle_history_arr = symbolic_ket_evolution_sp(M, l)
    println(L"The initial state $\ket{%$l,s} after %$M roundtrips has evolved into the following contributions:")
end