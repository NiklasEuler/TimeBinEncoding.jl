
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

    pops_fs_real, pops_fs_imag = pops_fs_phase_estimation(ρ_pure)
    relative_phases = initial_state_phase_estimation(pops_pure, pops_fs_real, pops_fs_imag)
    ρ_corrected = phase_on_density_matrix(ρ_pure, -1 * relative_phases)

    angles_compound_all = angles_compound(N)
    pops_fs_all_pure = fs_pop_compound(ρ_corrected, angles_compound_all)

    @test isapprox(
        coherence_extraction_compound(pops_pure, pops_fs_all_pure), 1.0, atol = 1e-8
    )

    ρ_nophase = density_matrix(
        insert_initial_state(correlated_timebin_state((2 + 3 * im) * ones(N)))
    )
    pops_nophase = populations(ρ_nophase)

    pops_fs_real, pops_fs_imag = pops_fs_phase_estimation(ρ_nophase)
    relative_phases =
        initial_state_phase_estimation(pops_nophase, pops_fs_real, pops_fs_imag)
    ρ_nocorrect = phase_on_density_matrix(ρ_nophase, -1 * relative_phases)

    @test relative_phases ≈ zeros(Float64, N)
    @test ρ_nocorrect ≈ ρ_nophase

    φ_arr = (0:N - 1) .* (π / 2)
    phase_gradient = cis.(φ_arr)
    ρ_phase_gradient = density_matrix(
        insert_initial_state(correlated_timebin_state(phase_gradient))
    )
    ρ_real = density_matrix(
        insert_initial_state(correlated_timebin_state(ones(N)))
    )

    pops_phase_gradient = populations(ρ_phase_gradient)

    pops_fs_real, pops_fs_imag = pops_fs_phase_estimation(ρ_phase_gradient)
    relative_phases =
        initial_state_phase_estimation(pops_phase_gradient, pops_fs_real, pops_fs_imag)
    ρ_correct = phase_on_density_matrix(ρ_phase_gradient, -1 * relative_phases)

    @test ρ_correct ≈ ρ_real

end


@testset "populations" begin
    N = 8
    ϵ = 0.1
    n_samples = 1e7

    wf_coeffs = cis.(2 * rand(N) * π)
	Ψ_init = insert_initial_state(correlated_timebin_state(wf_coeffs))
	ρ_mixed = density_matrix(Ψ_init)

    pops_mixed = populations(ρ_mixed)
    @test sum(pops_mixed) ≈ 1
    @test diag(ρ_mixed) ≈ pops_mixed
    @test typeof(pops_mixed) <: Vector{Float64}

    pops_mixed_sampled = populations(ρ_mixed, n_samples)
    @test isapprox(pops_mixed_sampled, pops_mixed, atol = 1e-3)
    @test sum(pops_mixed_sampled) ≈ 1

    @test isapprox(sample_populations(0.67, n_samples), 0.67, atol = 1e-3)

    pops_non_normalized = [0.1, 0.2, 0.3, 0.1]
    pops_non_normalized_sampled = sample_populations(
        pops_non_normalized, n_samples; unity = false
    )
    @test length(pops_non_normalized_sampled) == length(pops_non_normalized) == 4
    @test isapprox(pops_non_normalized_sampled, pops_non_normalized, atol = 1e-3)

    pops_to_much = [0.1, 0.8, 0.2]
    @test_throws ArgumentError sample_populations(pops_to_much, n_samples)

    pops_negative = [-0.1, 0.8, 0.2]
    @test_throws ArgumentError sample_populations(pops_negative, n_samples)


end
