using Test
using TimeBinEncoding
using LinearAlgebra

@testset "lcmk <-> j" begin
    # Write your tests here.
    N = 5
    j_arr = 1:(4*N^2)
    lcmk_arr = zeros(4, 4*N^2)
    j_arr_reconstructed = zeros(4*N^2)
    for j in eachindex(j_arr)
       lcmk = collect(j2lcmk(N,j))
       lcmk_arr[:,j] = lcmk
       j_arr_reconstructed[j] = lcmk2j(N,lcmk...)
    end
    @test j_arr == j_arr_reconstructed
    @test lcmk_arr[:,1] == [0,0,0,0]
    @test lcmk_arr[:,end] == [N-1,1,N-1,1]
    @test lcmk_arr[:,77] == [3,1,3,0]
    @test_throws ArgumentError j2lcmk(0,1)
    @test_throws ArgumentError j2lcmk(N,0)
    @test_throws ArgumentError j2lcmk(N,4*N^2+1)
    @test_throws ArgumentError lcmk2j(N,4,0,N,0)
    @test_throws ArgumentError lcmk2j(N,N,0,2,0)
    @test_throws ArgumentError lcmk2j(N,4,0,-1,0)
    @test_throws ArgumentError lcmk2j(N,-2,0,2,0)
    @test_throws ArgumentError lcmk2j(N,4,100,2,0)
    @test_throws ArgumentError lcmk2j(N,1,0,2,500)
    @test_throws ArgumentError lcmk2j(N,3,0,2,-1)
    @test_throws ArgumentError lcmk2j(N,2,-1,2,1)
    
end

@testset "lm <-> j" begin
    # Write your tests here.
    N = 5
    j_arr = 1:(N^2)
    lm_arr = zeros(2, N^2)
    j_arr_reconstructed = zeros(N^2)
    for j in eachindex(j_arr)
       lm = collect(j2lm(N,j))
       lm_arr[:,j] = lm
       j_arr_reconstructed[j] = lm2j(N,lm...)
    end
    @test j_arr == j_arr_reconstructed
    @test lm_arr[:,1] == [0,0]
    @test lm_arr[:,end] == [N-1,N-1]
    @test lm_arr[:,13] == [2,2]
    @test_throws ArgumentError j2lm(0,1)
    @test_throws ArgumentError j2lm(N,0)
    @test_throws ArgumentError j2lm(N,N^2+1)
    @test_throws ArgumentError lm2j(N,4,N)
    @test_throws ArgumentError lm2j(N,N,2)
    @test_throws ArgumentError lm2j(N,4,-1)
    @test_throws ArgumentError lm2j(N,-2,2)
    
end

@testset "lc <-> j" begin
    # Write your tests here.
    N = 5
    j_arr = 1:(2*N)
    lc_arr = zeros(2, 2*N)
    j_arr_reconstructed = zeros(2*N)
    for j in eachindex(j_arr)
       lc = collect(j2lc(j))
       lc_arr[:,j] = lc
       j_arr_reconstructed[j] = lc2j(lc...)
    end
    @test j_arr == j_arr_reconstructed
    @test lc_arr[:,1] == [0,0]
    @test lc_arr[:,end] == [N-1,1]
    @test lc_arr[:,6] == [2,1]
    @test_throws ArgumentError j2lc(0)
    @test_throws ArgumentError lc2j(-1,1)
    @test_throws ArgumentError lc2j(1,-1)
    @test_throws ArgumentError lc2j(4,2)
end

@testset "beam_splitter_operator" begin
    @test beam_splitter_operator(0) ≈ [[1,0] [0,1]]
    @test beam_splitter_operator(π/2) ≈ [[0,im] [im,0]]
    @test beam_splitter_operator(π/4) ≈ 1/√2 *[[1,im] [im,1]]
end

@testset "symbolic_final_state_projection_sp" begin

    @test_throws ArgumentError symbolic_final_state_projection_sp(0,0,0)
    @test_throws ArgumentError symbolic_final_state_projection_sp(0,0)
    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs = symbolic_final_state_projection_sp(1, 0, 0)
    @test j_idx_arr_fs == [1]
    @test trigonometric_history_arr_fs == [fill(0,1,1)]
    @test trigonometric_history_arr_fs == [fill(0,1,1)]
    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs = symbolic_final_state_projection_sp(2, 2, 1)
    @test j_idx_arr_fs == [1, 3]
    @test trigonometric_history_arr_fs == [[[1] [0]], [[0] [1]]]
    @test angle_history_arr_fs == [[[0] [1]], [[1] [1]]]
    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs = symbolic_final_state_projection_sp(2, 0, 1)
    @test j_idx_arr_fs == trigonometric_history_arr_fs == angle_history_arr_fs == []
end

@testset "explicit_state_evolution" begin
 	N = 2
 	M = 2
 	angles_1 = [0.5,0]*π
 	angles_2 = [0,0.25,0]*π
    angles = [angles_1,angles_2]

    ϕ = 0.0
	wf_coeffs = [cis(2*n*ϕ*π) for n in 0:N-1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)

    Ψ_out_numerical = mesh_evolution(Ψ_init, angles)
    Ψ_out_analytical = explicit_state_evolution(Ψ_init, angles)
    @test Ψ_out_numerical ≈ Ψ_out_analytical

    N = 4
	M = 4

	angles_1 = [0.25,0.25,0.25,0.25]*π
	angles_2 = [0.25,0.25,0.25,0.25,0.25]*π
	angles_3 = [0,0,0,0,0,0]*π
	angles_4 = [0.25,0.25,0.25,0.25,0.25,0.25,0.25]*π
    angles = [angles_1,angles_2,angles_3,angles_4]

    #ϕ = 0.3
	#wf_coeffs = [cis(2*n*ϕ*π) for n in 0:N-1]
    wf_coeffs = cis.(2 * rand(N) * π)
    tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)

    Ψ_out_numerical = mesh_evolution(Ψ_init, angles)
    Ψ_out_analytical = explicit_state_evolution(Ψ_init, angles)
    @test Ψ_out_numerical ≈ Ψ_out_analytical
end

@testset "explicit_final_state_projection_sp" begin
	angles_1 = [0.5,0]*π
	angles_2 = [0,0.25,0]*π

	# angles_1 = [0.25,0.25]*π
	# angles_2 = [0.25,0.25,0.25]*π
	
	angles = [angles_1,angles_2] 
    @test explicit_final_state_projection_sp(0,0, angles) == ([], [])
    @test explicit_final_state_projection_sp(0,1, angles) == ([], [])
    @test all(explicit_final_state_projection_sp(1,0, angles) .≈ ([1,3], 1/√2 .* Complex.([-1,1])))
    @test explicit_final_state_projection_sp(1,1, angles) == ([], [])
    @test explicit_final_state_projection_sp(2,0, angles) == ([], [])
    @test all(explicit_final_state_projection_sp(2,1, angles) .≈ ([1,3], im/√2 .* Complex.([1,1])))
    @test explicit_final_state_projection_sp(3,0, angles) == ([], [])
    

end

@testset "explicit_final_state_projection" begin
	N=2
    M=2
    
    # angles_1 = [0.5,0]*π
	# angles_2 = [0,0.25,0]*π

	angles_1 = [0.25,0.25]*π
	angles_2 = [0.25,0.25,0.25]*π
	
	angles = [angles_1,angles_2]
    j = lcmk2j(N+M,0,0,0,0)
    @test all(explicit_final_state_projection(j, angles) .≈ ([1], 1/4 .* Complex.([1])))
    j = lcmk2j(N+M,1,0,1,0)
    @test all(explicit_final_state_projection(j, angles) .≈ ([1,3,9,11], 1/4 .* Complex.([1,-1,-1,1])))
    j = lcmk2j(N+M,2,1,2,1)
    @test all(explicit_final_state_projection(j, angles) .≈ ([1,3,9,11], 1/4 .* Complex.([-1,-1,-1,-1])))
    j = lcmk2j(N+M,2,1,0,0)
    @test all(explicit_final_state_projection(j, angles) .≈ ([1,9], im/4 .* Complex.([1,1])))

end

@testset "density_matrix" begin
    N = 6
    ϕ = 0.64
	wf_coeffs = [cis(2*n*ϕ*π) for n in 0:N-1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
    
    ρ = density_matrix(Ψ_init)
    @test ρ == ρ'
    trace = diag(ρ)
    @test sum(trace) ≈ 1
    @test all(Float64.(trace) .≥ 0)
    Ψ_unnorm = 2(1+im) .* copy(Ψ_init)
    ρ_unnorm = density_matrix(Ψ_unnorm)
    trace = diag(ρ_unnorm)
    @test sum(trace) ≈ 1
    @test all(Float64.(trace) .≥ 0)
end


@testset "density_matrix_dephased" begin
    N = 6
    ϕ = 0.64
    ϵ = 0.1
	wf_coeffs = [cis(2*n*ϕ*π) for n in 0:N-1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
    
    ρ = density_matrix_dephased(Ψ_init,ϵ)
    @test ρ == ρ'
    trace = diag(ρ)
    @test sum(trace) ≈ 1
    @test all(Float64.(trace) .≥ 0)
end

@testset "phase_on_density_matrix" begin
    N = 6
    ϕ = 0.64
	wf_coeffs = [cis(2*n*ϕ*π) for n in 0:N-1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
    
    ρ = density_matrix(Ψ_init)
    φ_arr = zeros(N)
    ρ_rot = phase_on_density_matrix(ρ, φ_arr)
    @test ρ == ρ_rot
    φ_arr = 2*π*rand(N)
    ρ_rot = phase_on_density_matrix(ρ, φ_arr)
    @test ρ_rot == ρ_rot'
    ρ_rot_tor = phase_on_density_matrix(ρ_rot, -1*φ_arr)
    @test ρ_rot_tor ≈ ρ

end

@testset "phase_correction" begin
    N = 8
    wf_coeffs = cis.(2 * rand(N) * π)
	Ψ_init = insert_initial_state(correlated_timebin_state(wf_coeffs))
	ρ_pure = density_matrix(Ψ_init)
    ρ_corrected, relative_phases_auto = initial_state_phase_estimation(ρ_pure)
    @test isapprox(compound_coherence_extraction(ρ_corrected), 1.0, atol=1e-8)
    ρ_nophase = density_matrix(insert_initial_state(correlated_timebin_state((2+3im)*ones(N))))
	ρ_nocorrect, relative_phases = initial_state_phase_estimation(ρ_nophase)
    @test relative_phases ≈ zeros(Float64, N)
    @test ρ_nocorrect ≈ ρ_nophase
end

@testset "compound coherence extraction" begin
    N = 8
    ϵ = 0.1
    ϵ_angles = 0.05 * π
    wf_coeffs = cis.(2 * rand(N) * π)
	Ψ_init = insert_initial_state(correlated_timebin_state(wf_coeffs))
	ρ_pure = density_matrix(Ψ_init)
    Ψ_mes =  insert_initial_state(correlated_timebin_state(fill(1/sqrt(N),N)))
	mes_fidelity = fidelity(Ψ_mes, ρ_pure)
    @test isapprox(compound_coherence_extraction(ρ_pure), mes_fidelity, atol=1e-8)
    ρ_mixed = density_matrix_dephased(Ψ_init, ϵ)
    @test compound_coherence_extraction(ρ_mixed) ≤  compound_coherence_extraction(ρ_pure)
    @test compound_coherence_extraction(ρ_mixed, ϵ_angles) ≤ compound_coherence_extraction(ρ_mixed)
end
