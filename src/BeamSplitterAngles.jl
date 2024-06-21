export angles_kth_neighbor_interference, noisy_angles_symmetric
export angles_phase_estimation, angles_compound, angles_single_setup
export angles4bins_01, angles4bins_02, angles4bins_03, angles4bins
export angles_ref_bin_all_pairs, angles_pairs_from_mask
export graph_coloring


"""
    angles_kth_neighbor_interference(N, k)
    angles_kth_neighbor_interference(N, k, ϵ_angles)

Return an Vector of complete angle settings to exclusively interfere all combination of two
time bins with distance k.

Can also be called with optional argument `ϵ_angles`, which adds uniform random noise to the
angles.

See also [`noisy_angles_symmetric`](@ref).
"""
function angles_kth_neighbor_interference end

function angles_kth_neighbor_interference(N, k)
    N = convert(Int64, N)::Int64
    k = convert(Int64, k)::Int64

    @argcheck N ≥ 1
    @argcheck k ≥ 1
    @argcheck k < N

    angles = Vector{Vector{Vector{Float64}}}(undef, N - k)
    for i in 1:N - k
        angles[i] = _angles_kth_neighbor(N, k, i)
    end

    return angles
end

function angles_kth_neighbor_interference(N, k, ϵ_angles)
    N = convert(Int64, N)::Int64
    k = convert(Int64, k)::Int64

    @argcheck N ≥ 1
    @argcheck k ≥ 1
    @argcheck k < N

    angles = Vector{Vector{Vector{Float64}}}(undef, N - k)
    for i in 1:N - k
        angles[i] = _angles_kth_neighbor(N, k, i, ϵ_angles)
    end

    return angles
end

function _angles_kth_neighbor(N, k, i)
    angles_k_i = [zeros(Float64, n) for n in N:N + k]
    angles_k_i[1][i] = 0.5  * π # send the early time bin into the long loop
    angles_k_i[end][i+k] = 0.25 * π
    # interfere after k loops / k time bins travelled
    return angles_k_i
end

function _angles_kth_neighbor(N, k, i, ϵ_angles)
    angles_k_i = _angles_kth_neighbor(N, k, i)
    angles_k_i_noisy = noisy_angles_symmetric(angles_k_i, ϵ_angles)
    return angles_k_i_noisy
end


"""
    angles_phase_estimation(N::Real)
    angles_phase_estimation(N::Real, ϵ_angles)
    angles_phase_estimation(ρ_init::AbstractMatrix)

Return the set of all beam-splitter configurations `angles` needed for phase estimation.

Optionally, `ϵ_angles` can be given as an argument, returning uniform random angles
distributed symmetrically around the targeted angles. Futhermore, for convenience, the
function can instead of `N` also be called with a `ρ_init` argument to simply the use in
default argument usage.

# Returns

- `angles`::Vector{Vector{Vector{Float64}}}: Nested Vector of beam-splitter angles.
At the lowest two levels, it contains complete sets of beam splitter configurations as
Vector{Vector{Float64}}. All configurations that interfere initial-state time bins of
distance `k = 1` are collected into Vector{Vector{Vector{Float64}}}.

See also [`angles_single_setup`](@ref), [`angles_compound`](@ref), [`j_out_compound`](@ref),
[`coherence_extraction_compound`](@ref), [`noisy_angles_symmetric`](@ref).
"""
function angles_phase_estimation end

function angles_phase_estimation(N::Real)
    N = convert(Int64, N)::Int64
    angles = angles_kth_neighbor_interference(N, 1) # nearest neighbor phase estimation
    return angles
end

function angles_phase_estimation(N::Real, ϵ_angles)
    N = convert(Int64, N)::Int64
    angles = angles_kth_neighbor_interference(N, 1, ϵ_angles)
    # noisy nearest neighbor phase estimation
    return angles
end

function angles_phase_estimation(ρ_init::AbstractMatrix)
    # call signature to simplify pops_fs_phase_estimation default argument call signature
    N = ρ2N(ρ_init)
    angles =  angles_phase_estimation(N::Real)
    return angles
end

"""
    angles_compound(N)
    angles_compound(N, ϵ_angles)

Return the set of all beam-splitter configurations `angles` needed for compound coherence
extraction.

Optionally, `ϵ_angles` can be given as an argument, returning uniform random angles
distributed symmetrically around the targeted angles.

# Returns

- `angles`::Vector{Vector{Vector{Vector{Float64}}}}: Nested Vector of beam-splitter angles.
At the lowest two levels, it contains complete sets of beam splitter configurations as
Vector{Vector{Float64}}. All configurations that interfere initial-state time bins of
distance `k` are collected into Vector{Vector{Vector{Float64}}}, which are then grouped into
a final Vector{Vector{Vector{Vector{Float64}}}} for all allowed values of `k`

See also [`j_out_compound`](@ref), [`coherence_extraction_compound`](@ref),
[`noisy_angles_symmetric`](@ref), [`angles_single_setup`](@ref).
"""
function angles_compound end

function angles_compound(N)
    N = convert(Int64, N)::Int64

    angles = [angles_kth_neighbor_interference(N, k) for k in 1:N - 1]
    return angles:: Vector{Vector{Vector{Vector{Float64}}}}
end

function angles_compound(N, ϵ_angles)
    N = convert(Int64, N)::Int64
    ϵ_angles = convert(Float64, ϵ_angles)::Float64

    angles = [angles_kth_neighbor_interference(N, k, ϵ_angles) for k in 1:N - 1]
    return angles::Vector{Vector{Vector{Vector{Float64}}}}
end


"""
    noisy_angles_symmetric(angles, ϵ_angles)

Given a scheduled beam-splitter angle configuration `angles`, compute a noisy randomly
distributed actual realization of that configuration.

Each noisy angle `̂φ_i` is drawn symmetrically around its planned angle `φ_i`, i.e.,
`̂φ_i ∼ U([φ_i-ϵ_angles, φ_i+ϵ_angles])`.

See also [`initial_state_phase_estimation`](@ref).
"""
function noisy_angles_symmetric(angles, ϵ_angles)
    angles = convert(Vector{Vector{Float64}}, angles)::Vector{Vector{Float64}}
    ϵ_angles = convert(Float64, ϵ_angles)::Float64
    noisy_angles = [angles_m .+ ϵ_angles * 2 .* (rand(length(angles_m)) .- 1)
        for angles_m in angles]
    return noisy_angles
end

"""
    angles_single_setup(N)

Return the angles for a specific beam-splitter setup resulting in full coherence information
in every single readout time bin.

The setup requires `M = 2(N - 1)` round trips and interferes every input state time bin with
every other time bin in each of the non-vanishing output time bins.
"""
function angles_single_setup(N)
    N = convert(Int64, N)::Int64
    @argcheck isinteger(log2(N))
    M = 2 * (N - 1) # number of round trips
    N_half = N ÷ 2

    angles_cascade = [zeros(Float64, n) for n in N:N + M - 1]
    angles_cascade[1][1:N_half] .= π / 2 # send early half into long loop

    _recursive_beam_splitter_array!(
        N_half, 1 + N_half, 1 + N_half, angles_cascade, "center"
    )

    for i in 2:N_half - 1 # 2 are already included in minimum structure
        angles_cascade[N + i * 2][N + i] = π / 4 # interfere early and late branches
    end

    return angles_cascade
end

function _recursive_beam_splitter_array!(N_bs, n_idx, m_idx, angles, branch)
    @argcheck branch in ["early", "center", "late"]

    angles[m_idx][n_idx:n_idx + N_bs - 1] .= π / 4 # put original bs
    if N_bs == 1 # can happen in the flanks, no further recursion
        if branch == "early"
            angles[m_idx + 1][[n_idx]] .= π / 2 # left flank transparent couplers
        elseif branch == "late"
            angles[m_idx + 1][[n_idx + 1]] .= π / 2 # right flank transparent couplers
        end
    elseif N_bs == 2
        # smallest regular structure in the center, also appears in the flanks.
        # no further recursion
        angles[m_idx + 1][[n_idx, n_idx + 2]] .= π / 2
        angles[m_idx + 1][[n_idx + 1]] .= π / 4
        angles[m_idx + 3][[n_idx + 2]] .= π / 4
        if branch == "early"
            angles[m_idx + 4][[n_idx + 1]] .= π / 2 # left flank transparent couplers
            angles[m_idx + 4][[n_idx + 2]] .= π / 2 # left flank transparent couplers
        elseif branch == "late"
            angles[m_idx + 4][[n_idx + 3]] .= π / 2 # left flank transparent couplers
            angles[m_idx + 4][[n_idx + 4]] .= π / 2 # left flank transparent couplers
        end
    else
        N_bs_half = Int64(N_bs / 2)
        N_bs_quarter = Int64(N_bs / 4)
        angles[m_idx + N_bs_half][n_idx:n_idx + N_bs_quarter - 1] .= π / 2
        # left flank transparent couplers
        angles[m_idx + N_bs_half][
            n_idx + N_bs + N_bs_quarter:n_idx + N_bs + N_bs_half - 1] .= π / 2
        # right flank transparent couplers
        _recursive_beam_splitter_array!(
            N_bs_quarter, n_idx + N_bs_quarter, m_idx + N_bs_half + N_bs_quarter,
            angles, "early"
        ) # left flank beam splitter array
        _recursive_beam_splitter_array!(
            N_bs_quarter, n_idx + N_bs+N_bs_quarter, m_idx + N_bs_half + N_bs_quarter,
            angles, "late"
        ) # right flank beam splitter array
        _recursive_beam_splitter_array!(
            N_bs_half, n_idx + N_bs_half, m_idx + N_bs_half, angles, "center"
        ) # center branch beam splitter array
   end

    return nothing
end

"""
    angles4bins_01(N, l, m, p, q)

Compute the beam-splitter angles for interference of 4 bins using the first configuration,
in which the first and second time bins and the third and fourth time bins are interfered
with each other before the two pairs are interfered with each other.

# Arguments
- `N`: The number of time bins.
- `l`: The index of the first time bin.
- `m`: The index of the second time bin.
- `p`: The index of the third time bin.
- `q`: The index of the fourth time bin.

# Returns
An array of angles for interference of 4 bins using the first configuration.

See also [`angles4bins_02`](@ref), [`angles4bins_03`](@ref), [`angles4bins`](@ref).

"""
function angles4bins_01(N, l, m, p, q)
    angles = [zeros(i) for i in N:N + q - l]
    angles[1][[l + 1, p + 1]] .= 0.5
    angles[m - l + 1][m + 1] = 0.25
    angles[q - p + 1][q + 1] = 0.25
    angles[end][q + 1] = 0.25
    angles *= π
    return angles
end

"""
    angles4bins_02(N, l, m, p, q)

Compute the beam-splitter angles for interference of 4 bins using the second configuration,
in which the first and third time bins and the second and fourth time bins are interfered
with each other before the two pairs are interfered with each other.

# Arguments
- `N`: The number of time bins.
- `l`: The index of the first time bin.
- `m`: The index of the second time bin.
- `p`: The index of the third time bin.
- `q`: The index of the fourth time bin.

# Returns
An array of angles for interference of 4 bins using the second configuration.

See also [`angles4bins_01`](@ref), [`angles4bins_03`](@ref), [`angles4bins`](@ref).


"""
function angles4bins_02(N, l, m, p, q)
    angles = [zeros(i) for i in N:N + q - l]
    angles[1][[l + 1, m + 1]] .= 0.5
    angles[p - l + 1][p + 1] = 0.25
    angles[q - m + 1][q + 1] = 0.25
    angles[end][q + 1] = 0.25
    angles *= π
    return angles
end

"""
    angles4bins_03(N, l, m, p, q)

Compute the beam-splitter angles for interference of 4 bins using the  third configuration,
in which the first and fourth time bins and the second and third time bins are interfered
with each other before the two pairs are interfered with each other.

# Arguments
- `N`: The number of time bins.
- `l`: The index of the first time bin.
- `m`: The index of the second time bin.
- `p`: The index of the third time bin.
- `q`: The index of the fourth time bin.

# Returns
An array of angles for interference of 4 bins using the third configuration.

See also [`angles4bins_01`](@ref), [`angles4bins_02`](@ref), [`angles4bins`](@ref).


"""
function angles4bins_03(N, l, m, p, q)
    angles = [zeros(i) for i in N:N + q - l + 1]
    angles[1][l + 1] = 0.5
    angles[m - l + 2][m + 1] = 0.5
    angles[p - l + 2][p + 1] = 0.25
    angles[q - l + 1][q + 1] = 0.25
    angles[end][q + 1] = 0.25
    angles *= π
    return angles
end

"""
    angles4bins(N, l, m, p, q)

Compute the angles for time bin encoding with 4 bins using all three configurations.

# Arguments
- `N`: The number of time bins.
- `l`: The index of the first time bin.
- `m`: The index of the second time bin.
- `p`: The index of the third time bin.
- `q`: The index of the fourth time bin.

# Returns
A tuple of arrays of angles for interference of 4 bins using all three configurations.

See also [`angles4bins_01`](@ref), [`angles4bins_02`](@ref), [`angles4bins_03`](@ref).

"""
function angles4bins(N, l, m, p, q)
    angles_01 = angles4bins_01(N, l, m, p, q)
    angles_02 = angles4bins_02(N, l, m, p, q)
    angles_03 = angles4bins_03(N, l, m, p, q)

    return angles_01, angles_02, angles_03
end

function angles_ref_bin_all_pairs(N, idx_ref; population_bins=false)
	N = convert(Int64, N)::Int64 # initial number of time bins
    idx_ref = convert(Int64, idx_ref)::Int64 # reference time-bin index
	@argcheck idx_ref >= 0
	@argcheck idx_ref < N

	Δfirst = idx_ref
	Δlast = N - 1 - idx_ref

	n_branch_vert = max(Δfirst - 1, 0)
        # number of secondary vertical branches going of the main diagonal branch
	n_branch_diag = max(Δlast - 1, 0)
        # number of secondary diagonal branches going of the main vertical branch

	M = max(Δfirst + n_branch_diag, Δlast + n_branch_vert) + 1 # depth of the mesh lattice

	angles = [zeros(Float64, n) for n in N:N + M - 1]
	offset_last = Δlast > 0 ? 1 : 0
        # correction term needed when there are bins later than the reference bin
	offset_first = Δfirst > 0 ? 1 : 0
        # correction term needed when there are bins earlier than the reference bin

	angles[1][idx_ref + 1] = asin(sqrt((n_branch_vert + offset_last)/(Δfirst + Δlast))) / π
	    # need enough time to split of enough packages for interference
        # with all bins later then the reference
	angles[n_branch_diag + 1][1:idx_ref] .= population_bins ? θ_pop_ref : 0.5
        # insert all bins earlier then the reference bin into the long loop
        # if population_bins is true, the beam is split to generate a population reference
        # output for further measurement options
	if n_branch_vert > 0
		angles[1][idx_ref + 2:end] .= 0.5
            # delay later time bins to make space for interference with earlier time bins
		angles[1 + n_branch_vert][idx_ref + 2 + n_branch_vert:end] .=
            population_bins ? θ_pop_ref : 0.5
            # reinsert into short bin. if population_bins is true, the beam is split to
            # generate a population reference output for further measurement options
		for n in 1:n_branch_vert
            # lso include potential mirror at the end of main diagonal branch
 			angles[1 + n][idx_ref + 1 + n] =
                asin(sqrt(1/(n_branch_vert + 1 + offset_last - n))) / π
			    # branch of ever larger amounts for overall homogeneous distro.
                # last entry is just 0.5
		end
	elseif population_bins
        angles[1][idx_ref + 2:end] .= θ_pop_ref
            # if population_bins is true, the beam is split to generate a population
            # reference output for further measurement options
    end

	for n in 1:n_branch_diag
        # also include potential mirror at the end of main vertical branch
		angles[1 + n][idx_ref + 1] =
            asin(sqrt(1/(n_branch_diag + 1 + offset_first - n))) / π
		    # branch off ever larger amplitudes for overall homogeneous distro.
            # last entry is just 0.5
	end

	idx_start_final_bs = Δfirst > 0 ? idx_ref + 1 : 2
        # idx where the earliest and the reference bin interfere
	angles[end][idx_start_final_bs:idx_start_final_bs + Δfirst + Δlast - 1] .= 0.25
        # actual time-bin interference
	angles .*= π

	return angles
end

function angles_pairs_from_mask(N, pair_arr; population_bins=false)
	N = convert(Int64, N)::Int64
	sort!.(pair_arr)
	incoming_diag_idxs = [pair[2] + 1 for pair in pair_arr]
	Δpairs = [pair[2] - pair[1] for pair in pair_arr]
	Δmax = max(Δpairs...)
	M = Δmax + 1
	angles = [zeros(Float64, n) for n in N:N + M - 1]
	roundtrip_insertion_idxs = Δmax .- Δpairs .+ 1
	for (pair_idx, pair) in enumerate(pair_arr)
		idx_early = pair[1] + 1
		idx_late = pair[2] + 1
		angles[roundtrip_insertion_idxs[pair_idx]][idx_early] =
            population_bins ? θ_pop_ref : 0.5
		angles[end][idx_late] = 0.25
	end
	if population_bins
		pair_arr_inv = sort(pair_arr, by=last, rev=true)
		for pair in pair_arr_inv
			idx_late = pair[2] + 1
			flag_population_bin = false
			for pop_branch_delay in 0:Δmax - 1
				pop_bin_idx = idx_late + Δmax - pop_branch_delay
				if !(pop_bin_idx in incoming_diag_idxs)
					angles[1 + pop_branch_delay][idx_late] = θ_pop_ref
					push!(incoming_diag_idxs, pop_bin_idx)
					flag_population_bin = true
					break
				end
			end
			if !flag_population_bin
                throw(ArgumentError(
                    "A population bin could not be established due to a lack of a timeslot."
                    )
                )
            end
        end
	end
    angles .*= π
	return angles
end

"""
    graph_coloring(N)

Given an integer `N`, this function performs graph-edge coloring to generate pairings of
bins/nodes. If `N` is even, it calls the `_graph_coloring_even` function to generate the
pairings. If `N` is odd, it first calls the `graph_coloring` function recursively with
`N + 1` to generate pairings for an even number of bins/nodes. Then, it removes the first
pair from each perfect matching to leave out one bin/node, resulting in an odd number of
bins/nodes.

# Arguments
- `N::Int`: The number of bins/nodes.

# Returns
- `pairings::Vector{Vector{Vector{Int}}}`: A vector of perfect matchings of the fully
    connected graph of `N` nodes. Each perfect matching is a vector containing the chosen
    edges of the matching each represented by the a vector with the two adjacent node
    indices of the respective edge.

"""
function graph_coloring(N)
    if iseven(N)
        pairings = _graph_coloring_even(N)
    else
        pairings_even = graph_coloring(N + 1)
        pairings = [perfect_matching[2:end] for perfect_matching in pairings_even]
        # delete the first pair, which is always the reference bin. This way, one bin is
        # left out, and the number of bins/nodes is odd again.
    end

    return pairings
end


function _graph_coloring_even(N)
	N = convert(Int64, N)::Int64
	@argcheck iseven(N)
	@argcheck N > 1

	N_half = N ÷ 2
	n_pairs = N_half - 1 # number of pairings beside the defining pair adjacent to node N-1
	pairings_all = [[[0, 0] for _ in 1:N_half] for _ in 1:N - 1]

	for i in 1:N - 1
		pairings_all[i][1][1] = i - 1 # first pair determined by pair between
		pairings_all[i][1][2] = N - 1 # reference and i'th bin
		for j in 1:n_pairs
			# symmetrically pairing up nodes with ever increasing distance to the defining
            # node. The boundaries of the chain are taken to be open boundary conditions
			val_left = mod(i - 1 - j, N - 1)
                # index of the node to the left of the defining node
			val_right = mod(i - 1 + j, N - 1)
                # index of the node to the right of the defining node
			pairings_all[i][1 + j] = sort([val_left, val_right])
		end
	end

	return pairings_all
end
