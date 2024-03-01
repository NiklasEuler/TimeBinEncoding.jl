using Test
using TimeBinEncoding
using LinearAlgebra
using SparseArrays

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

@testset "explicit_ket_evolution_sp" begin
    N = 2
 	M = 2
 	angles_1 = [0.5,0]*π
 	angles_2 = [0,0.25,0]*π
    angles = [angles_1,angles_2]

    wf_coeffs_norm = normalize([0,1])
	Ψ_sp = insert_initial_state_sp(wf_coeffs_norm)

    vec_mesh = mesh_evolution_sp(Ψ_sp, angles)
    idx, val = explicit_ket_evolution_sp(1, angles)
    vec_explicit = sparsevec(idx, val, 2*(N+M))
    @test vec_explicit ≈ vec_mesh

    N = 4
	M = 4

	angles_1 = [0.25,0.25,0.25,0.25]*π
	angles_2 = [0.25,0.25,0.25,0.25,0.25]*π
	angles_3 = [0,0,0,0,0,0]*π
	angles_4 = [0.25,0.25,0.25,0.25,0.25,0.25,0.25]*π
    angles = [angles_1,angles_2,angles_3,angles_4]

    wf_coeffs_norm = normalize([0,0,1,0])
	Ψ_sp = insert_initial_state_sp(wf_coeffs_norm)

    vec_mesh = mesh_evolution_sp(Ψ_sp, angles)
    idx, val = explicit_ket_evolution_sp(2, angles)
    vec_explicit = sparsevec(idx, val, 2*(N+M))
    @test vec_explicit ≈ vec_mesh
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

    ψ_dense = zeros(ComplexF64, 4*N^2)
	ψ_dense[1] = 1
	ψ_sparse = spzeros(ComplexF64, 4*N^2)
	ψ_sparse[1] = 1
    @test mesh_evolution(ψ_dense, angles) ≈ mesh_evolution(ψ_sparse, angles) 

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


    N = 4
    M = 2*(N-1)
    angles = angles_single_setup(N)
    for j in 1:TimeBinEncoding.n_loops2*(N+M)^2
        j_idx_arr_contr_symb, coeff_arr_symb = TimeBinEncoding.explicit_final_state_projection_symbolic_backend(N, M, j, angles)
        j_idx_arr_contr_mesh, coeff_arr_mesh = TimeBinEncoding.explicit_final_state_projection_mesh(N, M, j, angles)
        @test j_idx_arr_contr_symb == j_idx_arr_contr_mesh
        @test coeff_arr_symb ≈ coeff_arr_mesh
    end
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

@testset "purity" begin
    N = 6
    ϕ = 0.64
    ϵ = 0.1

	wf_coeffs = [cis(2*n*ϕ*π) for n in 0:N-1]
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
    @test not_pure_at_all ≈ 1/N^2

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

@testset "angles_single_setup" begin
    N = 4
    angles_1_mod = [0.5,0.5,0,0]*π
	angles_2_mod = [0,0,0,0,0]*π
	angles_3_mod = [0,0,0.25,0.25,0,0]*π
	angles_4_mod = [0,0,0.5,0.25,0.5,0,0]*π
	angles_5_mod = [0,0,0,0,0,0,0,0]*π
	angles_6_mod = [0,0,0,0,0.25,0,0,0,0]*π
	angles_mod = [angles_1_mod,angles_2_mod,angles_3_mod,angles_4_mod,angles_5_mod,angles_6_mod]
    angles = angles_single_setup(N)
    @test angles == angles_mod

    N = 8
    M_mod = 14
	angles_8_mod = [zeros(Float64, n) for n in N:N+M_mod-1]
	angles_8_mod[1][1:4] .= 0.5*π
	angles_8_mod[5][5:8] .= 0.25*π
	angles_8_mod[7][[5,10]] .= 0.5*π
	angles_8_mod[7][7:8] .= 0.25*π
	angles_8_mod[8][[6,8,10]] .= 0.25*π
	angles_8_mod[8][[7,9]] .= 0.5*π
	angles_8_mod[9][[6,11]] .= 0.5*π
	angles_8_mod[10][[9]] .= 0.25*π
	angles_8_mod[12][10] = 0.25*π
	angles_8_mod[14][11] = 0.25*π
    angles_8 = angles_single_setup(N)
    @test angles_8 == angles_8_mod

    N = 16
    M_mod = 30
    angles_16_mod_symm = [zeros(Float64, n) for n in N:N+M_mod-1]
	angles_16_mod_symm[1][1:8] .= 0.5*π
	angles_16_mod_symm[9][9:16] .= 0.25*π
	angles_16_mod_symm[13][[9,10,19,20]] .= 0.5*π
	angles_16_mod_symm[13][13:16] .= 0.25*π
	angles_16_mod_symm[15][[13,18]] .= 0.5*π
	angles_16_mod_symm[15][[11,12,15,16,19,20]] .= 0.25*π
	angles_16_mod_symm[16][[11,13,15,17,19,21]] .= 0.5*π
	angles_16_mod_symm[16][[12,14,16,18,20]] .= 0.25*π
	angles_16_mod_symm[17][[14,19]] .= 0.5*π
	angles_16_mod_symm[18][[13,17,21]] .= 0.25*π
	angles_16_mod_symm[19][[12,13,22,23]] .= 0.5*π
	angles_16_mod_symm[20][[18]] .= 0.25*π
	angles_16_mod_symm[22][[19]] .= 0.25*π
	angles_16_mod_symm[24][[20]] .= 0.25*π
	angles_16_mod_symm[26][[21]] .= 0.25*π
	angles_16_mod_symm[28][[22]] .= 0.25*π
	angles_16_mod_symm[30][[23]] .= 0.25*π
    angles_16 = angles_single_setup(N)
    @test angles_16 == angles_16_mod_symm

    N = 6
    @test_throws ArgumentError angles_single_setup(N)
end

@testset "j_out_single_setup" begin
    N = 8
    M_mod = 14
    j_mod_arr = [lcmk2j(N+M_mod,7,0,7,0), lcmk2j(N+M_mod,14,1,14,1), lcmk2j(N+M_mod,8,0,8,0), lcmk2j(N+M_mod,13,1,13,1),
    lcmk2j(N+M_mod,9,0,9,0), lcmk2j(N+M_mod,12,1,12,1), lcmk2j(N+M_mod,10,0,10,0), lcmk2j(N+M_mod,11,1,11,1)]
    j_arr = j_out_single_setup(N)
    @test all([j ∈ j_arr for j in j_mod_arr])
    @test all([j ∈ j_mod_arr for j in j_arr])
    N = 10
    @test_throws ArgumentError j_out_single_setup(N)
end