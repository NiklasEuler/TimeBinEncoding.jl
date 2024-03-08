
@testset "density_matrix" begin
    N = 6
    ϕ = 0.64
	wf_coeffs = [cis(2 * n * ϕ * π) for n in 0:N - 1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
    ρ = density_matrix(Ψ_init)
    @test ρ == ρ'
    trace = diag(ρ)
    @test sum(trace) ≈ 1
    @test all(Float64.(trace) .≥ 0)
    Ψ_unnorm = 2 * (1 + im) .* copy(Ψ_init)
    ρ_unnorm = density_matrix(Ψ_unnorm)
    trace = diag(ρ_unnorm)
    @test sum(trace) ≈ 1
    @test all(Float64.(trace) .≥ 0)
end


@testset "density_matrix_dephased" begin
    N = 6
    ϕ = 0.64
    ϵ = 0.1
	wf_coeffs = [cis(2 * n * ϕ * π) for n in 0:N - 1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
    ρ = density_matrix_dephased(Ψ_init,ϵ)
    @test ρ == ρ'
    trace = diag(ρ)
    @test sum(trace) ≈ 1
    @test all(Float64.(trace) .≥ 0)
end

@testset "purity" begin
    N = 6
    ϕ = 0.64
    ϵ = 0.1
	wf_coeffs = [cis(2 * n * ϕ * π) for n in 0:N - 1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
    ρ = density_matrix(Ψ_init)

    pure = purity(ρ)
    @test pure ≈ 1
    @test typeof(pure) == Float64

    ρ_mixed = density_matrix_dephased(Ψ_init, ϵ)
    not_so_pure = purity(ρ_mixed)
    @test not_so_pure < 1.0

    ρ_inf_temp = density_matrix_dephased(Ψ_init, 1)
    not_pure_at_all = purity(ρ_inf_temp)
    @test not_pure_at_all ≈ 1 / N^2

end

@testset "phase_on_density_matrix" begin
    N = 6
    ϕ = 0.64
	wf_coeffs = [cis(2 * n * ϕ * π) for n in 0:N - 1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
    ρ = density_matrix(Ψ_init)

    φ_arr = zeros(N)
    ρ_rot = phase_on_density_matrix(ρ, φ_arr)
    @test ρ == ρ_rot

    φ_arr = 2 * π*rand(N)
    ρ_rot = phase_on_density_matrix(ρ, φ_arr)
    @test ρ_rot == ρ_rot'

    ρ_rot_tor = phase_on_density_matrix(ρ_rot, -1 * φ_arr)
    @test ρ_rot_tor ≈ ρ

end

@testset "phase_correction" begin
    N = 8
    wf_coeffs = cis.(2 * rand(N) * π)
	Ψ_init = insert_initial_state(correlated_timebin_state(wf_coeffs))
	ρ_pure = density_matrix(Ψ_init)
    pops_pure = populations(ρ_pure)

    ρ_corrected, relative_phases_auto = initial_state_phase_estimation(ρ_pure)

    angles_compound_all = compound_angles(N)
    pops_fs_all_pure = pops_fs_compound(ρ_corrected, angles_compound_all)

    @test isapprox(
        compound_coherence_extraction(pops_pure, pops_fs_all_pure), 1.0, atol = 1e-8
    )
    ρ_nophase = density_matrix(
        insert_initial_state(correlated_timebin_state((2 + 3 * im) * ones(N)))
    )
	ρ_nocorrect, relative_phases = initial_state_phase_estimation(ρ_nophase)
    @test relative_phases ≈ zeros(Float64, N)
    @test ρ_nocorrect ≈ ρ_nophase
end


@testset "populations" begin
    N = 8
    ϵ = 0.1
    n_samples = 1e6

    wf_coeffs = cis.(2 * rand(N) * π)
	Ψ_init = insert_initial_state(correlated_timebin_state(wf_coeffs))
	ρ_mixed = density_matrix(Ψ_init)

    pops_mixed = populations(ρ_mixed)
    @test sum(pops_mixed) ≈ 1
    @test diag(ρ_mixed) .≈ pops_mixed
    @test typeof(pops_mixed) <: Vector{Float64}
    pops_mixed_sampled = populations(ρ_mixed, n_samples)

    @test isapprox(pops_mixed_sampled, pops_mixed, atol = 1e-3)
end
