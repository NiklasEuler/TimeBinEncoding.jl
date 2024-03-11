export symbolic_ket_evolution_sp
export symbolic_fs_projection_sp

"""
    symbolic_ket_evolution_sp(M, l)
    symbolic_ket_evolution_sp(M)

Symbolically compute the evolution of the single-photon ket `|l,S⟩` through `M` roundtrips
in the `|l,c⟩` basis.

If the argument `l` is omitted, the evolution of `|0,S⟩` is returned.

# Returns

 -`j_idx_arr`::Vector{Int64}: All `j` indices of the `|l,c⟩` for which the state after
    evolution has a potentially non-zero coefficient.
 -`trigonometric_history_arr`::Vector{Matrix{Int64}}: Contains a Matrix for each index `j`
    from `j_idx_arr` which keeps track at what times cos or sin function should be applied.
    Each of the columns of each Matrix corresponds to one round trip in the fiber-loop
    system. Each row marks the contributions of a single time-bin trajectory through the
    lattice. A entry of 0 indicates a cos factor for that given round-trip index and
    trajectory, whereas a value of 1 indicates a sin contribution, i.e., the photon changes
    from the short to the long loop or vice versa.
 -`angle_history_arr`:Vector{Matrix{Int64}}: Of the same shape as
    `trigonometric_history_arr`. Every entry holds information about the temporal position
    of each trajectory at each round trip. Each value corresponds to the current time-bin
    index `l` in the `|l,c⟩` basis at each round trip.

    See also [`symbolic_fs_projection_sp`](@ref), [`explicit_ket_evolution_sp`](@ref),
    [`explicit_ket_evolution`](@ref).

"""
function symbolic_ket_evolution_sp end

function symbolic_ket_evolution_sp(M, l)
    M = convert(Int64, M)::Int64 # number of roundtrips
    l = convert(Int64, l)::Int64 # initial time bin index
    #@argcheck M > 0
    @argcheck l ≥ 0

    j_idx_arr, trigonometric_history_arr, angle_history_arr = symbolic_ket_evolution_sp(M)
    # compute evolution of |0,S⟩ to be adapted for
    # |l,S⟩ by adding l to all time-bin indices.
    if l ≠ 0
        for i in eachindex(angle_history_arr)
            angle_history_arr[i] .+= l
        end
        j_idx_arr .+= N_LOOPS * l
   end

    return j_idx_arr, trigonometric_history_arr, angle_history_arr
end

function symbolic_ket_evolution_sp(M)
    M = convert(Int64, M)::Int64 # number of roundtrips

    @argcheck M > 0

    trigonometric_history_arr = fill(Matrix{Int64}(undef, 0, M), N_LOOPS * (M + 1))
    # initiate at full length
    trigonometric_history_arr[1] = fill(-1, 1, M)
    # initiate with illegal value to track potential bugs
    angle_history_arr = copy(trigonometric_history_arr)
    j_idx_arr = [1, collect(1 + N_LOOPS:N_LOOPS * M)..., N_LOOPS * (M + 1)]
    for m in 1:M # current round trip index
        _symbolic_ket_coin_sp!(m, trigonometric_history_arr, angle_history_arr)
        # perform symbolic beam splitter application
        _symbolic_ket_shift_sp!(M, m, trigonometric_history_arr, angle_history_arr)
        # perform symbolic time-bin shift
    end
    deleteat!(trigonometric_history_arr, [N_LOOPS, length(trigonometric_history_arr) - 1])
    # remove undef matrices for the two empty bins
    deleteat!(angle_history_arr, [N_LOOPS, length(angle_history_arr) - 1])
    # remove undef matrices for the two empty bins

    return j_idx_arr, trigonometric_history_arr, angle_history_arr
end

function _symbolic_ket_coin_sp!(m, trigonometric_history_arr, angle_history_arr)
    m = convert(Int64, m)::Int64 # current round trip index

    @argcheck m > 0

    for j in 1:N_LOOPS:m*N_LOOPS
        l, c  = j2lc(j)
        trigonometric_history_short = copy(trigonometric_history_arr[j])
        # copy of trigonometric history of short loop bin to facilitate coupling
        trigonometric_history_long = copy(trigonometric_history_arr[j + 1])
        # copy of trigonometric history of long loop bin to facilitate coupling
        trigonometric_history_arr[j][:, m] .= 0 # cos to stay in short loop
        trigonometric_history_arr[j + 1][:, m] .= 0 # cos to stay in long loop
        trigonometric_history_short[:, m] .= 1 # sin to switch to long loop
        trigonometric_history_long[:, m] .= 1 # sin to switch to short loop
        trigonometric_history_arr[j] =
            vcat(trigonometric_history_arr[j], trigonometric_history_long)
        trigonometric_history_arr[j + 1] =
            vcat(trigonometric_history_arr[j + 1], trigonometric_history_short)

        angle_history_arr[j][:, m] .= l# l'th time bin in short loop
        angle_history_arr[j + 1][:, m] .= l # l'th time bin in short loop
        angle_history_short = copy(angle_history_arr[j])
        # copy of angle history of short loop bin to facilitate coupling
        angle_history_long = copy(angle_history_arr[j + 1])
        # copy of angle history of long loop bin to facilitate coupling
        angle_history_arr[j] = vcat(angle_history_arr[j], angle_history_long)
        # add new possible trajectories switching from long to short
        angle_history_arr[j + 1] = vcat(angle_history_arr[j + 1], angle_history_short)
        # add new possible trajectories switching from short to long
    end

    return nothing
end

function _symbolic_ket_shift_sp!(M, m, trigonometric_history_arr, angle_history_arr)
    M = convert(Int64, M)::Int64 # number of roundtrips
    m = convert(Int64, m)::Int64 # current round trip index

    @argcheck m > 0
    @argcheck m ≤ M
    @argcheck M > 0

    for j in N_LOOPS * (M + 1):-N_LOOPS:2 * N_LOOPS
        #go through vector from the back so that no values are overwritten when shifted
        trigonometric_history_arr[j] = trigonometric_history_arr[j - N_LOOPS]
        # shift cos/sin information
        angle_history_arr[j] = angle_history_arr[j - N_LOOPS]
        # shift time-bin trajectories
    end

    trigonometric_history_arr[N_LOOPS] = Matrix{Int64}(undef, 0, M)
    # reset 0L bin after shift
    angle_history_arr[N_LOOPS] = Matrix{Int64}(undef, 0, M)
    # reset 0L bin after shift
    return nothing
end


"""
    symbolic_fs_projection_sp(M, l, c)
    symbolic_fs_projection_sp(M, j)

Symbolically compute the projection of the single-photon ket `|l,c⟩` onto the final state in
terms of the initial state coefficients.

Accepts either the time-bin index `l` and loop index `c` or the `j` index of the final state
`|l,c⟩` in the `|l,c⟩` basis.

# Returns

-`j_idx_arr`::Vector{Int64}: All `j` indices of kets in the initial state which have nonzero
    contributions to `|l,c⟩` after evolution with the generic mesh lattice.
 -`trigonometric_history_arr`::Vector{Matrix{Int64}}: Contains a Matrix for each index `j`
    from `j_idx_arr` which keeps track at what times cos or sin function should be applied.
    Each of the columns of each Matrix corresponds to one round trip in the fiber-loop
    system. Each row marks the contributions of a single time-bin trajectory through the
    lattice. A entry of 0 indicates a cos factor for that given round-trip index and
    trajectory, whereas a value of 1 indicates a sin contribution, i.e., the photon changes
    from the short to the long loop or vice versa.
 -`angle_history_arr`:Vector{Matrix{Int64}}: Of the same shape as
    `trigonometric_history_arr`. Every entry holds information about the temporal position
    of each trajectory at each round trip. Each value corresponds to the current time-bin
    index `l` in the `|l,c⟩` basis at each round trip.

    See also [`symbolic_ket_evolution_sp`](@ref), [`explicit_fs_projection_sp`](@ref),
    [`explicit_fs_projection`](@ref).
"""
function symbolic_fs_projection_sp end

function symbolic_fs_projection_sp(M, l, c)
    M = convert(Int64, M)::Int64 # number of roundtrips
    l = convert(Int64, l)::Int64 # final state time bin index
    c = convert(Int64, c)::Int64 # final state loop index

    @argcheck M > 0

    return _symbolic_fs_projection_worker_sp(M, l, c)
end

function symbolic_fs_projection_sp(M, j)
    M = convert(Int64, M)::Int64 # number of roundtrips
    j = convert(Int64, j)::Int64 # final state ket index

    @argcheck M > 0

    l, c = j2lc(j)
    return _symbolic_fs_projection_worker_sp(M, l, c)
end

function _symbolic_fs_projection_worker_sp(M, l, c)
    trigonometric_history_arr_fs = Vector{Matrix{Int64}}(undef, 0)
    angle_history_arr_fs = copy(trigonometric_history_arr_fs)
    j_idx_arr_fs = Int64[]

    j_fs = lc2j(l, c )
    if j_fs ≠ 2 # no contribution to |0L> bin
        if c == 1
            startidx = max(l - M, 0) # can at most travel M time bins, or start at zero
            endidx = max(l - 1, 0)
            # cannot start at same index, as at least one switch to long loop in necessary
        else
            startidx = max(l - M + 1, 0)
            # has to switch back once, so can travel at most M-1 time bins, or start at zero
            endidx = l
            # can start at same index, as no switch to long loop is necessary
        end

        for l_init in startidx:endidx
            j_idx_arr, trigonometric_history_arr, angle_history_arr =
                symbolic_ket_evolution_sp(M, l_init) # all possible initial state kets
            idx = searchsortedfirst(j_idx_arr, j_fs)
            j_init = lc2j(l_init, 0) # only short loop initial states
            push!(j_idx_arr_fs, j_init)
            push!(trigonometric_history_arr_fs, trigonometric_history_arr[idx])
            push!(angle_history_arr_fs, angle_history_arr[idx])
        end
   end

    return j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs
end
