export correlated_timebin_state, insert_initial_state, insert_initial_state_sp
export density_matrix, density_matrix_dephased, white_noise
export fidelity, purity, populations, sample_populations
#export density_matrix_dephased_identical, white_noise_identical


"""
    correlated_timebin_state(wf_coeffs::Vector)

Convert `wf_coeffs` holding the coefficients of correlated two-photon time bin populations
to the `|l, m ⟩' basis.

See also [`insert_initial_state`](@ref), [`density_matrix`](@ref).
"""
function correlated_timebin_state(wf_coeffs::AbstractVector)
    wf_coeffs = convert(Vector{ComplexF64}, wf_coeffs)::Vector{ComplexF64}
    N = length(wf_coeffs)
    coeffs = normalize(wf_coeffs)
    time_bin_state_vec = zeros(ComplexF64, N^2)
    for i in 0:N - 1
        j = lm2j(N, i, i) # correlated |ii⟩ contributions only
        time_bin_state_vec[j] = coeffs[i + 1]
   end

    return time_bin_state_vec
end


"""
    insert_initial_state(time_bin_state_vec::Vector)

Insert the two-photon time-bin state in the `|l, m ⟩' basis into the short loop by imbedding
`|l, m ⟩' → `|l, 0, m , 0⟩'.

See also [`insert_initial_state_sp`](@ref), [`correlated_timebin_state`](@ref),
[`density_matrix`](@ref).
"""
function insert_initial_state(time_bin_state_vec::Vector)
    N = Int64(sqrt(length(time_bin_state_vec)))
    full_state_vec = zeros(ComplexF64, length(time_bin_state_vec) * N_LOOPS2)

    for l in 0:N - 1, m in 0:N - 1
        j_time_bin = lm2j(N, l, m)
        j_full = lcmk2j(N, l, 0, m, 0) #insertion into the short-short ket |l0m0⟩
        full_state_vec[j_full] = time_bin_state_vec[j_time_bin]
   end

    return full_state_vec
end

"""
    insert_initial_state_sp(time_bin_state_vec::Vector)

Insert the single-photon time-bin state in the `|l⟩' basis into the short loop by imbedding
`|l⟩' → `|l, 0⟩'.

See also [`insert_initial_state`](@ref), [`correlated_timebin_state`](@ref),
[`density_matrix`](@ref).
"""
function insert_initial_state_sp(time_bin_state_vec::Vector)
    N = Int64(length(time_bin_state_vec))
    full_state_vec = zeros(ComplexF64, length(time_bin_state_vec) * N_LOOPS)

    for l in 0:N - 1
        j_full = lc2j(l, 0) #insertion into the short ket |l0⟩
        full_state_vec[j_full] = time_bin_state_vec[l + 1]
   end

    return full_state_vec
end

"""
    density_matrix(Ψ)

Compute the density matrix `ρ` to the wave function `Ψ`.

See also `density_matrix_dephased`.
"""
function density_matrix(Ψ)
    Ψ = convert(Vector{ComplexF64}, Ψ)::Vector{ComplexF64}
    normalize!(Ψ) # normalization

    ρ = kron(Ψ, Ψ')
    return ρ
end

"""
    density_matrix_dephased(Ψ, ϵ)

Compute a dephased density matrix `ρ` to the wave function `Ψ`.

The dephasing is introduced as white noise in the short-short `|i0j0⟩` populations:
ρ = (1-ϵ) * |Ψ⟩⟨Ψ| + ϵ/N^2 * ∑_i, j |i0j0⟩⟨i0j0|

See also [`density_matrix`](@ref), [`white_noise`](@ref).
"""
function density_matrix_dephased(Ψ, ϵ)
    Ψ = convert(Vector{ComplexF64}, Ψ)::Vector{ComplexF64}
    ϵ = convert(Float64, ϵ)::Float64
    N = Int64(sqrt(length(Ψ)) / N_LOOPS)::Int64

    @argcheck ϵ ≥ 0
    @argcheck ϵ ≤ 1

    ρ_pure = density_matrix(Ψ)
    ρ = (1 - ϵ) * ρ_pure + ϵ * white_noise(N)
    return ρ
end

#= function density_matrix_dephased_identical(Ψ, ϵ)
    Ψ = convert(Vector{ComplexF64}, Ψ)::Vector{ComplexF64}
    ϵ = convert(Float64, ϵ)::Float64
    d_hilbert_space = Int(sqrt(length(Ψ)))

    N = Int(-1 / 4 + sqrt(1 / 16 + d_hilbert_space / 2)) # p-q formula

    @argcheck ϵ ≥ 0
    @argcheck ϵ ≤ 1

    ρ_pure = density_matrix(Ψ)
    ρ = (1 - ϵ) * ρ_pure + ϵ * white_noise_identical(N)
    return ρ
end =#

"""
    phase_on_density_matrix(ρ, φ_arr)

Apply phases `φ_arr` to the correlated time bins of the density matrix `ρ`

Returns a new density matrix `ρ` after phase application. Each time bin |ii⟩, i ∈ {0, 1,…}
is subjected to phase `φ_arr[i + 1]`.

See also [`initial_state_phase_estimation`](@ref).
"""
function phase_on_density_matrix(ρ, φ_arr)
    ρ_rot = convert(Matrix{ComplexF64}, copy(ρ))::Matrix{ComplexF64}
    φ_arr = convert(Vector{Float64}, φ_arr)::Vector{Float64}
    N = ρ2N(ρ)

    @argcheck length(φ_arr) == N # number of phases should match number of time bins

    tb_idxs = [lcmk2j(N, i, 0, i, 0) for i in 0:N - 1]
    # all valid correlated initial-state time bins
    for (idx1, j1) in enumerate(tb_idxs), (idx2, j2) in enumerate(tb_idxs)
        ρ_rot[j1, j2] *= cis(φ_arr[idx1] - φ_arr[idx2])
        # apply relative phase to every coherence
   end

    return ρ_rot
end

"""
    white_noise(N)

Compute a normalized dephased density matrix with equal populations for all short-short time
bins.

ρ_wn = 1 / N^2 * ∑_i, j |i0j0⟩⟨i0j0|

See also [`density_matrix_dephased`](@ref).
"""
function white_noise(N)
    ρ_noise = zeros(ComplexF64, N_LOOPS2 * N^2, N_LOOPS2 * N^2)
    weight = 1 / N^2
    for l in 0:N - 1, m in 0:N - 1
        j = lcmk2j(N, l, 0, m , 0)
        ρ_noise[j, j] = weight
   end

    return ρ_noise
end

"""
    fidelity(Ψ::Vector,ρ::Matrix)

Compute the fidelity between density matrix `ρ` and pure state `Ψ`.
"""
function fidelity(Ψ::AbstractVector,ρ::AbstractMatrix)
    Ψ = convert(Vector{ComplexF64}, Ψ)::Vector{ComplexF64}
    fidelity = Ψ' * ρ * Ψ
    return convert(Float64, real(fidelity))
end

"""
    purity(ρ)

Compute the purity of density matrix `ρ`.
"""
function purity(ρ)
    return Float64(real(sum(diag(ρ * ρ))))
end

"""
    ρ2N(ρ)

Extract the time bin number `N` from the density matrix `ρ`.
"""
function ρ2N(ρ)
    return Int64(sqrt(size(ρ)[1] / (N_LOOPS2)))
end

"""
    populations(ρ::Matrix)
    populations(ρ::Matrix, n_samples)

Return the state populations of the density matrix `ρ`.

Optionally, the populations can be subjected to finite statistics. This is achieved through
the optional `n_samples` function argument, which specifies the number of simulated samples
used for the reconstruction of the population distribution.

See also [`sample_populations`](@ref)
"""
function populations end

function populations(ρ::Matrix)
    pops = convert(Vector{Float64}, diag(ρ))::Vector{Float64}

    @argcheck sum(pops) ≈ 1 # unity check

    return pops
end

function populations(ρ::Matrix, n_samples)
    pops = populations(ρ::Matrix)
    measured_pops = sample_populations(pops, n_samples)
    return measured_pops
end

"""
    sample_populations(pop::Real, n_samples)
    sample_populations(pops::Vector{<:Real}, n_samples; unity=true)

Return an approximate normalized reconstruction of the original distribution `pops` based on
the measurement statistics generated with `n_samples` random samples.

Accepts also a single (unnormalized) argument `pop` as a standalone probability to sample
from. If `unity=true`, the sum of the probabilities is required to be unity, otherwise also
not normalized distributions are accepted.

See also [`populations`](@ref).
"""
function sample_populations end

function sample_populations(pop::Real, n_samples)
    pops = [pop, 1 - pop]
    pops_sampled = sample_populations(pops, n_samples)
    pop_sampled = pops_sampled[1]
    return pop_sampled
end

function sample_populations(pops::Vector{<:Real}, n_samples; unity=true)
    n_samples = convert(Int64, n_samples)::Int64 # number of measurements

    if(unity)
        @argcheck sum(pops) ≈ 1 # unity check
    else
        @argcheck 0 < sum(pops) < 1 # normalizable check
    end
    @argcheck n_samples > 0
    @argcheck all(pops .≥ 0)

    if !unity
        push!(pops, 1 - sum(pops)) # append the remaining probability to the last entry
    end

    pops_measured = zero(pops)
    samples_ordered =
        sample(collect(1:length(pops)), ProbabilityWeights(pops), n_samples, ordered=true)
        # vector of n_sample entries, each entry being the j index of the measured outcome
    count_map_samples = collect(countmap(samples_ordered))
    # vector of pairs, the first entry being the j index and the second the count number in
    # the simulated measurement

    for count_pair in count_map_samples
        pops_measured[count_pair.first] = count_pair.second/n_samples # normalization
    end

    if !unity
        pop!(pops) # remove the appended entry
        pop!(pops_measured)
    end

    return pops_measured
end

if !@isdefined cispi
    function cispi(x)
        # define the complex exponential function with π factor for older Julia versions
        return cis(x .* π)
    end
end
