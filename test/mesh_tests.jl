
@testset "beam_splitter_operator" begin
    @test beam_splitter_operator(0) ≈ [[1, 0] [0, 1]]
    @test beam_splitter_operator(π / 2) ≈ [[0, im] [im, 0]]
    @test beam_splitter_operator(π / 4) ≈ 1 / √2 *[[1, im] [im, 1]]
end

@testset "mesh_evolution" begin

    N = 2
 	M = 2
 	angles_1 = [0.5, 0] * π
 	angles_2 = [0, 0.25, 0] * π
    angles = [angles_1, angles_2]

    ϕ = 0
	wf_coeffs = [cis(2 * n * ϕ * π) for n in 0:N - 1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
    ρ_init = density_matrix(Ψ_init)

    Ψ_out_vector = mesh_evolution(Ψ_init, angles)
    ρ_out_vector = density_matrix(Ψ_out_vector)
    ρ_out_matrix = mesh_evolution(ρ_init, angles)
    @test ρ_out_vector ≈ ρ_out_matrix

    N = 8

    angles = angles_single_setup(N)

    #ϕ = 0.3
	#wf_coeffs = [cis(2 * n * ϕ * π) for n in 0:N - 1]
    wf_coeffs = cis.(2 * rand(N) * π)
    tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
    ρ_init = density_matrix(Ψ_init)

    Ψ_out_vector = mesh_evolution(Ψ_init, angles)
    ρ_out_vector = density_matrix(Ψ_out_vector)
    ρ_out_matrix = mesh_evolution(ρ_init, angles)
    @test ρ_out_vector ≈ ρ_out_matrix
end

@testset "mesh_evolution_sp" begin

    N = 2
 	M = 2
 	angles_1 = [0.1, 0.3] * π
 	angles_2 = [0.7, 0.25, 0.2] * π
    angles = [angles_1, angles_2]

    wf_coeffs_norm = normalize([1 + 0.2im, 0.5 + 0.6im])
	Ψ_init = insert_initial_state_sp(wf_coeffs_norm)
    ρ_init = density_matrix(Ψ_init)

    Ψ_out_vector = mesh_evolution_sp(Ψ_init, angles)
    ρ_out_vector = density_matrix(Ψ_out_vector)
    ρ_out_matrix = mesh_evolution_sp(ρ_init, angles)
    @test ρ_out_vector ≈ ρ_out_matrix

    N = 4

    angles = angles_single_setup(N)

    #ϕ = 0.3
	#wf_coeffs = [cis(2 * n * ϕ * π) for n in 0:N - 1]
    wf_coeffs_norm = normalize(cis.(2 * rand(N) * π))
	Ψ_init = insert_initial_state_sp(wf_coeffs_norm)
    ρ_init = density_matrix(Ψ_init)

    Ψ_out_vector = mesh_evolution_sp(Ψ_init, angles)
    ρ_out_vector = density_matrix(Ψ_out_vector)
    ρ_out_matrix = mesh_evolution_sp(ρ_init, angles)
    @test ρ_out_vector ≈ ρ_out_matrix
end
