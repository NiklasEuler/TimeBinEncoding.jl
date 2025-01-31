@testset "mesh_evolution_identical" begin
    N = 4

    d_local_hs_bl = N * (2 * N + 1)
    d_full_hs_bl = d_local_hs_bl ^ 2
        # hilbert-space dimension of 2 species a 2 photons in full short/long notion

    Ψ_init = spzeros(ComplexF64, d_full_hs_bl)
	j = lcmk2j_super_identical(N, 0, 0, 1, 0, 0, 0, 1, 0)
	Ψ_init[j] = 1

    angles = [zeros(n) for n in N:N+1]
	angles[1][1] = 0.5π
	angles[2][2] = 0.25π

    Ψ_out =  mesh_evolution_identical(Ψ_init, angles)
    @test sum(abs2, Ψ_out) ≈ 1 # normalization

    j1 = lcmk2j_super_identical(N + 2, 1, 0, 1, 0, 1, 0, 1, 0)
    j2 = lcmk2j_super_identical(N + 2, 1, 0, 1, 0, 2, 1, 2, 1)
    j3 = lcmk2j_super_identical(N + 2, 2, 1, 2, 1, 1, 0, 1, 0)
    j4 = lcmk2j_super_identical(N + 2, 2, 1, 2, 1, 2, 1, 2, 1)

    @test Ψ_out[[j1, j2, j3, j4]] ≈ [0.5, 0.5, 0.5, 0.5] # Hong-Ou-Mandel effect

    kron_mem = TimeBinEncoding._kron_mem_arr(N, 2)
    Ψ_out_mem =  mesh_evolution_identical(Ψ_init, angles, kron_mem)
    @test Ψ_out_mem ≈ Ψ_out

    Ψ_in_dense = Vector(Ψ_init)
    Ψ_out_dense = mesh_evolution_identical(Ψ_in_dense, angles)
    @test Ψ_out ≈ Ψ_out_dense

    ρ_out_manual = density_matrix(Ψ_out)
    ρ_in_manual = density_matrix(Ψ_init)
    ρ_out = mesh_evolution_identical(ρ_in_manual, angles)
    @test ρ_out ≈ ρ_out_manual

    ρ_out_mem = mesh_evolution_identical(ρ_in_manual, angles, kron_mem)
    @test ρ_out_mem ≈ ρ_out

    N = 5
    d_local_hs_bl = N * (2 * N + 1)

    Ψ_init_sp = spzeros(ComplexF64, d_local_hs_bl)
    j = lcmk2j_identical(N, 1, 0, 3, 0)
	Ψ_init_sp[j] = 1

    angles = [zeros(n) for n in N:N+2]
	angles[1][2] = 0.5π
	angles[3][4] = 0.25π

    Ψ_out_sp =  mesh_evolution_sp_identical(Ψ_init_sp, angles)
    @test sum(abs2, Ψ_out_sp) ≈ 1 # normalization

    j1 = lcmk2j_identical(N + 3, 3, 0, 3, 0)
    j2 = lcmk2j_identical(N + 3, 4, 1, 4, 1)

    @test Ψ_out_sp[[j1, j2]] ≈ [-1/sqrt(2), -1/sqrt(2)] # Hong-Ou-Mandel effect

    Ψ_in_dense_sp = Vector(Ψ_init_sp)
    Ψ_out_dense_sp = mesh_evolution_sp_identical(Ψ_in_dense_sp, angles)
    @test Ψ_out_sp ≈ Ψ_out_dense_sp

    ρ_out_manual_sp = density_matrix(Ψ_out_sp)
    ρ_in_manual_sp = density_matrix(Ψ_init_sp)
    ρ_out_sp = mesh_evolution_sp_identical(ρ_in_manual_sp, angles)
    @test ρ_out_sp ≈ ρ_out_manual_sp

    d_full_hs_bl = d_local_hs_bl ^ 2
    # hilbert-space dimension of 2 species a 2 photons in full short/long notion

    Ψ_init = spzeros(ComplexF64, d_full_hs_bl)
    j = lcmk2j_super_identical(N, 1, 0, 3, 0, 1, 0, 3, 0)
	Ψ_init[j] = 1

    Ψ_out =  mesh_evolution_identical(Ψ_init, angles)
    Ψ_out_sp_kron = kron(Ψ_out_sp, Ψ_out_sp)
    @test Ψ_out ≈ Ψ_out_sp_kron

    angles = [zeros(n) for n in N:N+3]

    Ψ_out_sp =  mesh_evolution_sp_identical(Ψ_init_sp, angles)

    N = 9
    d_local_hs_bl = N * (2 * N + 1)
    Ψ_out_sp_manual = spzeros(ComplexF64, d_local_hs_bl)
    j = lcmk2j_identical(N, 1, 0, 3, 0)
    Ψ_out_sp_manual[j] = 1
    @test Ψ_out_sp ≈ Ψ_out_sp_manual


end
