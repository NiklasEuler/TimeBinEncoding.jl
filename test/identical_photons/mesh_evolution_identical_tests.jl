@testset "mesh_evolution_identical" begin
    N = 4

    d_local_hs_b = Int(N * (N + 1) / 2) # hilbert-space dimension of two identical photons in N time bins. No notion of short/long
    d_local_hs_bl = N * (2 * N + 1)
    d_full_hs_b = d_local_hs_b ^ 2
    d_full_hs_bl = d_local_hs_bl ^ 2 # hilbert-space dimension of 2 species a 2 photons in full short/long notion

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

end
