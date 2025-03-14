export visualize_symbolic_ket_evolution_sp, visualize_symbolic_fs_projection_sp,
    visual_meas_coh_map, visualize_combined_measurement_coherence_map
    #export visualize_measurement_coherence_map_identical

"""
    visualize_symbolic_ket_evolution_sp(M, l_init)

Print the generic evolution of the single-photon initial state `|l_init, S⟩` after `M`
roundtrips. The function prints the contributing states and their weights computed
symbolically for general beam-splitter angles.

# Arguments
- `M`: The number of roundtrips.
- `l_init`: The initial state time-bin index.

# Returns
- `nothing`

See also [`symbolic_ket_evolution_sp`](@ref), ['visualize_symbolic_fs_projection_sp`](@ref).

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
        _trigonometric_string_formatter(trigonometric_history_arr[i], angle_history_arr[i])
        println("]")
   end

    return nothing
end




"""
    visualize_symbolic_fs_projection_sp(M, l_fs, c_fs)

Print the contributing single-photon initial state terms to a generic final state after `M`
roundtrips and subsequent projection onto `|l_fs, c_fs⟩⟨l_fs, c_fs|`. The function prints
the contributing states and their weights computed symbolically for general beam-splitter
angles.

# Arguments
- `M`: The number of roundtrips.
- `l_fs`: The value of `l_fs` in the final-state projector `|l_fs, c_fs⟩⟨l_fs, c_fs|`.
- `c_fs`: The value of `c_fs` in the final-state projector `|l_fs, c_fs⟩⟨l_fs, c_fs|`.

# Returns
- `nothing`

See also [`symbolic_fs_projection_sp`](@ref), ['visualize_symbolic_ket_evolution_sp`](@ref).

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
        _trigonometric_string_formatter(trigonometric_history_arr_fs[i],
            angle_history_arr_fs[i])
        println("]")
    end
    println("}")
    return nothing
end

function _trigonometric_string_formatter(trigonometric_history, angle_history)
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


"""
    visual_meas_coh_map(
        j_out::Int64,
        angles,
        phases=ones(Float64, length(angles[1]))::Vector,
        off_l=0,
        off_m=0;
        extract_diagonal=true
    )
    visual_meas_coh_map(
        j_out::Vector{Int64},
        angles,
        projector_weights=ones(Float64, length(j_out)),
        phases=ones(Float64, length(angles[1]))::Vector,
        off_l=0,
        off_m=0;
        extract_diagonal=true
    )

Visualize the contributing initial-state coherences contributing to the measured projectors
defined through `j_out`. The function prints the contributing coherences and their weights.
Optionally, initial state `phases` can be provided. The function also allows for time shifts
for the signal and idler photons through the `off_l` and `off_m` arguments.

# Arguments
- `j_out`: The output time bin(s) of the measured projector(s).
- `angles`: The beam-splitter angles along the interference measurement.
- `phases`: The phases of the roundtrips. Default is an array of ones.
- `off_l::Int64`: The time shift for the signal photon. Default is 0.
- `off_m::Int64`: The time shift for the signal photon. Default is 0.
- `extract_diagonal::Bool`: Whether to extract the diagonal elements. Default is true.

# Returns
- `nothing`

"""
function visual_meas_coh_map end

function visual_meas_coh_map(
    j_out::Int64,
    angles,
    phases=ones(Float64, length(angles[1]))::AbstractVector,
    off_l=0,
    off_m=0;
    extract_diagonal=true
)
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    projector_weight = 1 # needed for single projector for interoperability
    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map(j_out, angles, projector_weight, phases)

    _print_header_offsets(off_l, off_m)

    l_out, c_out, m_out, k_out = j2lcmk(N + M, j_out)
	println("⟨", l_out, " ", c_out, " ", m_out, " ", k_out, "|(SC)^M ρ (C^†S^†)^M) |",
        l_out, " ", c_out, " ", m_out, " ", k_out, "⟩ ="
    )
    _visual_coh(N, j1_arr, j2_arr, weights, extract_diagonal, phases, off_l, off_m)
    return nothing
end

function visual_meas_coh_map(
    j_out_arr::AbstractVector{Int64},
    angles,
    projector_weights=ones(Float64, length(j_out_arr)),
    phases=ones(Float64, length(angles[1]))::AbstractVector,
    off_l=0,
    off_m=0;
    extract_diagonal=true,
)
    M = length(angles)  # number of roundtrips
    N = length(angles[1]) # initial number of time bins

    j1_arr, j2_arr, weights =
        explicit_fs_coherence_map(j_out_arr, angles, projector_weights, phases)

    _print_header_offsets(off_l, off_m)

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
    _visual_coh(N, j1_arr, j2_arr, weights, extract_diagonal, phases, off_l, off_m)

    return nothing
end

function _print_header_offsets(off_l, off_m) # print header line with offset information
    if off_l != 0 || off_m != 0
        println("with offsets: off_l = ", off_l, ", off_m = ", off_m)
    end

    return nothing
end

function _visual_coh(
    N,
    j1_arr,
    j2_arr,
    weights,
    extract_diagonal=true,
    phases=ones(Float64, N),
    off_l=0,
    off_m=0
) # encoding of coherence information and printout
    contr_j_idxs = correlated_short_bins_idxs(N)
    extractable_correlated_coherences = []
    if phases == ones(Float64, N)
        display_weights = round.(Real.(weights), digits=5)
    else
        display_weights = round.(weights, digits=5)
    end


	for i in eachindex(j1_arr)
        j1 = j1_arr[i]
        j2 = j2_arr[i]
		l1, c1, m1, k1 = j2lcmk(N, j1)
		l2, c2, m2, k2 = j2lcmk(N, j2)

        l1 -= off_l
        l2 -= off_l
        m1 -= off_m
        m2 -= off_m

        if l1 < 0 || l2 < 0 || m1 < 0 || m2 < 0
            continue
        end
        if l1 >= N - off_m || l2 >= N - off_m
            continue
        end
        if m1 >= N - off_l || m2 >= N - off_l
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
