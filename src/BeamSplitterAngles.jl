export angles_kth_neighbor_interference, noisy_angles_symmetric
export angles_phase_estimation, angles_compound, angles_single_setup
export angles4bins_01, angles4bins_02, angles4bins_03, angles4bins
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
