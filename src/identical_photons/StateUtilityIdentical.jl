export density_matrix_dephased_identical, white_noise_identical
export krauss_operator_single_species, krauss_operators_dual_species
export photon_loss_channel, density_matrix_dephased_krauss_identical

function density_matrix_dephased_identical(Ψ, ϵ)
    Ψ = convert(Vector{ComplexF64}, Ψ)::Vector{ComplexF64}
    ϵ = convert(Float64, ϵ)::Float64
    d_hilbert_space = Int(sqrt(length(Ψ)))

    N = Int(-1 / 4 + sqrt(1 / 16 + d_hilbert_space / 2)) # p-q formular

    @argcheck ϵ ≥ 0
    @argcheck ϵ ≤ 1

    ρ_pure = density_matrix(Ψ)
    ρ = (1 - ϵ) * ρ_pure + ϵ * white_noise_identical(N)
    return ρ
end


"""
    white_noise_identical(N)

Generate a white noise density matrix for a four photon system with two pairs of identical
photons. This model allows for white noise populations where one of the photons is detected
in a wrong bin, which can differ by one time bin compared to the correct bin.

"""
function white_noise_identical(N)
    d_local_hs_bl = N * (2 * N + 1)
    # local Hilbert space dimension of two identical photons in the bin-loop basis
    ρ_diag = zeros(ComplexF64, d_local_hs_bl^2)
    for l in 0:N - 1, m in l:N - 1
        j = lcmk2j_super_identical(N, l, 0, m, 0, l, 0, m, 0)
        ρ_diag[j] = 1

        if l != m && l < N - 1
            j = lcmk2j_super_identical(N, l + 1, 0, m, 0, l, 0, m, 0)
            ρ_diag[j] = 1
            j = lcmk2j_super_identical(N, l, 0, m, 0, l + 1, 0, m, 0)
            ρ_diag[j] = 1
        end

        if m < N - 1
            j = lcmk2j_super_identical(N, l, 0, m + 1, 0, l, 0, m, 0)
            ρ_diag[j] = 1
            j = lcmk2j_super_identical(N, l, 0, m, 0, l, 0, m + 1, 0)
            ρ_diag[j] = 1
        end
   end

   ρ_diag /= sum(ρ_diag)
   ρ_noise = Diagonal(ρ_diag)

   return ρ_noise
end

function krauss_operator_single_species(N, l, m)
    N = convert(Int64, N)::Int64
    l = convert(Int64, l)::Int64
    m = convert(Int64, m)::Int64

    @argcheck N > 0
    @argcheck 0 ≤ l < N
    @argcheck 0 ≤ m < N

    d_local_hs_bl = N * (2 * N + 1)
    # local Hilbert space dimension of two identical photons in the bin-loop basis
    krauss_op = spzeros(Float64, d_local_hs_bl, d_local_hs_bl)
    for n in 0:N - 1
        bin1_init, bin2_init = sort([n, m])
        j1 = lcmk2j_identical(N, bin1_init, 0, bin2_init, 0)
        bin1_after, bin2_after = sort([n, l])
        j2 = lcmk2j_identical(N, bin1_after, 0, bin2_after, 0)
        if l ≠ n ≠ m
            coeff = 1
        elseif n == l == m
            coeff = 2
        else
            coeff = √2
        end
        krauss_op[j2, j1] = coeff
    end
    krauss_op *= 1 / (2 * √(N + 1))

    return krauss_op
end

function krauss_operators_dual_species(N, l, m)
    d_local_hs_bl = N * (2 * N + 1)
    # local Hilbert space dimension of two identical photons in the bin-loop basis
    I_local = Diagonal(ones(d_local_hs_bl))
    krauss_op = krauss_operator_single_species(N, l, m)
    klm_signal = kron(krauss_op, I_local)
    klm_idler = kron(I_local, krauss_op)

    return klm_signal, klm_idler
end

function photon_loss_channel(N, ρ)
    d_local_hs_bl = N * (2 * N + 1)
    d_full_hs_bl  = d_local_hs_bl ^ 2
    noisy_state = spzeros(ComplexF64, d_full_hs_bl, d_full_hs_bl)
    for l in 0:N - 1, m in 0:N - 1
        klm_signal, klm_idler = krauss_operators_dual_species(N, l, m)
        noisy_state += klm_signal * ρ * klm_signal' + klm_idler * ρ * klm_idler'
    end
    return noisy_state

end

function density_matrix_dephased_krauss_identical(N, ρ, ϵ)
    ρ_noisy = photon_loss_channel(N, ρ)
    ρ_dephased = (1 - ϵ) * ρ + ϵ * ρ_noisy
    return ρ_dephased
end

#= function white_noise_dark_count_identical(N, loss_factor)
    d_hilbert_space = N_LOOPS * N * (N_LOOPS * N + 1) / 2 # correct?
    ρ_diag = zeros(ComplexF64, d_hilbert_space^2)
    for l in 0:N - 1
        for m in l:N - 1, for k in l:N - 1
            j = lcmk2j_super_identical(N, l, 0, m, 0, k, 0, l, 0)
            ρ_diag[j] = 1
        end

end
 =#


#= function white_noise_identical(N)
    d_hilbert_space = N * (2 * N + 1)
    #ρ_noise = zeros(ComplexF64, d_hilbert_space^2, d_hilbert_space^2)
    ρ_diag = zeros(ComplexF64, d_hilbert_space^2)
    for l in 0:N - 1, m in l:N - 1
        for p in 0:N - 1, q in p:N - 1
            j = lcmk2j_super_identical(N, l, 0, m, 0, p, 0, q, 0)
            ρ_diag[j] = 1
        end
   end
   ρ_diag /= sum(ρ_diag)
   ρ_noise = Diagonal(ρ_diag)

   return ρ_noise
end =#
