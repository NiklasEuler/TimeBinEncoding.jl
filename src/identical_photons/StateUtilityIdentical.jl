export density_matrix_dephased_identical, white_noise_identical

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
    d_hilbert_space = N * (2 * N + 1)
    # simplified expression for N_LOOPS = 2 canceling out geometric series
    #ρ_noise = zeros(ComplexF64, d_hilbert_space^2, d_hilbert_space^2)
    ρ_diag = zeros(ComplexF64, d_hilbert_space^2)
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
