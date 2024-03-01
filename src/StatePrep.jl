export correlated_timebin_state, insert_initial_state, insert_initial_state_sp
export density_matrix, density_matrix_dephased, white_noise
export fidelity, purity



"""
    correlated_timebin_state(wf_coeffs::Vector)

Convert `wf_coeffs` holding the coefficients of correlated two-photon time bin popluations to the `|l,m⟩' basis.

See also `insert_initial_state`, `density_matrix`.
"""
function correlated_timebin_state(wf_coeffs::Vector)
    wf_coeffs = convert(Vector{ComplexF64}, wf_coeffs)::Vector{ComplexF64}
    N = length(wf_coeffs)
    coeffs = normalize(wf_coeffs)
    time_bin_state_vec = zeros(ComplexF64, N^2)
    for i in 0:N-1
        j = lm2j(N,i,i) # correlated |ii⟩ contributions only
        time_bin_state_vec[j] = coeffs[i+1]
    end
    return time_bin_state_vec
end


"""
    insert_initial_state(time_bin_state_vec::Vector)

Insert the two-photon time-bin state in the `|l,m⟩' basis into the short loop by imbedding `|l,m⟩' → `|l,0,m,0⟩'.

See also `insert_initial_state_sp`, `correlated_timebin_state`, `density_matrix`.
"""
function insert_initial_state(time_bin_state_vec::Vector)
    N = Int64(sqrt(length(time_bin_state_vec)))
    full_state_vec = zeros(ComplexF64, length(time_bin_state_vec)*n_loops2)
    for l in 0:N-1, m in 0:N-1
        j_time_bin = lm2j(N, l, m)
        j_full = lcmk2j(N, l, 0, m, 0) #insertion into the short-short ket |l0m0⟩
        full_state_vec[j_full] = time_bin_state_vec[j_time_bin]
    end
    return full_state_vec
end

"""
    insert_initial_state_sp(time_bin_state_vec::Vector)

Insert the single-photon time-bin state in the `|l⟩' basis into the short loop by imbedding `|l⟩' → `|l,0⟩'.

See also `insert_initial_state` `correlated_timebin_state`, `density_matrix`.
"""
function insert_initial_state_sp(time_bin_state_vec::Vector)
    N = Int64(length(time_bin_state_vec))
    full_state_vec = zeros(ComplexF64, length(time_bin_state_vec)*n_loops)
    for l in 0:N-1
        j_full = lc2j(l, 0) #insertion into the short ket |l0⟩
        full_state_vec[j_full] = time_bin_state_vec[l+1]
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
ρ = (1-ϵ) * |Ψ⟩⟨Ψ| + ϵ/N^2 * ∑_i,j |i0j0⟩⟨i0j0|

See also `density_matrix`, `white_noise`.
"""
function density_matrix_dephased(Ψ, ϵ)
    Ψ = convert(Vector{ComplexF64}, Ψ)::Vector{ComplexF64}
    ϵ = convert(Float64, ϵ)::Float64
    N = Int64(sqrt(length(Ψ))/n_loops)::Int64
    @argcheck ϵ ≥ 0
    @argcheck ϵ ≤ 1
    
    ρ_pure = density_matrix(Ψ)
    ρ = (1-ϵ) * ρ_pure + ϵ * white_noise(N) 
    return ρ
end


"""
    white_noise(N)

Compute a normalized dephased density matrix with equal populations for all short-short time bins.

ρ_wn = 1/N^2 * ∑_i,j |i0j0⟩⟨i0j0|
See also `density_matrix_dephased`.
"""
function white_noise(N)
    ρ_noise = zeros(ComplexF64, n_loops2*N^2,n_loops2*N^2)
    weight = 1/N^2
    for l in 0:N-1, m in 0:N-1
        j = lcmk2j(N,l,0,m,0)
        ρ_noise[j,j] = weight
    end
    return ρ_noise
end

"""
    fidelity(Ψ::Vector,ρ::Matrix)

Compute the fidelity between density matrix `ρ` and pure state `Ψ`.
"""
function fidelity(Ψ::Vector,ρ::Matrix)
    Ψ = convert(Vector{ComplexF64}, Ψ)::Vector{ComplexF64}
    fidelity = Ψ' * ρ * Ψ
    return convert(Float64,real(fidelity))
end

"""
    purity(ρ)

Compute the purity of density matrix `Ψ`
"""
function purity(ρ)
    return Float64(real(sum(diag(ρ*ρ))))
end

