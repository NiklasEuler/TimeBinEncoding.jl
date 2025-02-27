@testset "density_matrix_dephased_identical" begin
    N = 4

    d_local_hs_bl = N * (2 * N + 1)
    d_full_hs_bl = d_local_hs_bl ^ 2
        # hilbert-space dimension of 2 species a 2 photons in full short/long notion

    Random.seed!(4441)

    Ψ_init = spzeros(ComplexF64, d_full_hs_bl)
    phases = cis.(2 * rand(N) * π)
    for l in 0:N - 1
		for m in l:N - 1
			j_super = lcmk2j_super_identical(N, l, 0, m, 0, l, 0, m, 0)
			Ψ_init[j_super] = phases[l + 1] * phases[m + 1]
		end
	end
	normalize!(Ψ_init)
    ρ_pure = density_matrix(Ψ_init)

    ϵ = 0
    ρ_dephased = density_matrix_dephased_identical(Ψ_init, ϵ)
    @test ρ_dephased ≈ ρ_pure

    ρ_krauss = density_matrix_dephased_krauss_identical(N, ρ_pure, ϵ)
    @test ρ_krauss ≈ ρ_dephased

    ϵ = 0.1
    ρ_dephased = density_matrix_dephased_identical(Ψ_init, ϵ)
    @test ρ_dephased ≈ (1 - ϵ) * ρ_pure + ϵ * white_noise_identical(N)
    fidelity_deph = fidelity(Ψ_init, ρ_dephased)
    @test isapprox(fidelity_deph, 0.9017241379310346,atol=1e-10)
    # old implementation: 0.90294117647

    ρ_krauss = density_matrix_dephased_krauss_identical(N, ρ_pure, ϵ)
    fidelity_krauss = fidelity(Ψ_init, ρ_krauss)
    @test isapprox(fidelity_krauss, 0.91, atol=1e-10)
    @test fidelity_krauss ≥ fidelity_deph

    @test sum(diag(ρ_krauss)) ≈ 1
    @test sum(diag(ρ_dephased)) ≈ 1

end
