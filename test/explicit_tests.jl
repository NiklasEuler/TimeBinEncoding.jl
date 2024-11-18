@testset "explicit_ket_evolution_sp" begin
    N = 2
 	M = 2
 	angles_1 = [0.5, 0] * π
 	angles_2 = [0, 0.25, 0] * π
    angles = [angles_1, angles_2]

    wf_coeffs_norm = normalize([0, 1])
	Ψ_sp = insert_initial_state_sp(wf_coeffs_norm)

    vec_mesh = mesh_evolution_sp(Ψ_sp, angles)
    idx, val = explicit_ket_evolution_sp(1, angles)
    vec_explicit = sparsevec(idx, val, 2 * (N + M))
    @test vec_explicit ≈ vec_mesh

    N = 4
	M = 4

	angles_1 = [0.25, 0.25, 0.25, 0.25] * π
	angles_2 = [0.25, 0.25, 0.25, 0.25, 0.25] * π
	angles_3 = [0, 0, 0, 0, 0, 0] * π
	angles_4 = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25] * π
    angles = [angles_1, angles_2, angles_3, angles_4]

    wf_coeffs_norm = normalize([0, 0, 1, 0])
	Ψ_sp = insert_initial_state_sp(wf_coeffs_norm)

    vec_mesh = mesh_evolution_sp(Ψ_sp, angles)
    idx, val = explicit_ket_evolution_sp(2, angles)
    vec_explicit = sparsevec(idx, val, 2 * (N + M))
    @test vec_explicit ≈ vec_mesh
    vec_mesh_sparse = mesh_evolution_sp(sparse(Ψ_sp), angles)
    @test vec_explicit ≈ vec_mesh_sparse

end



@testset "explicit_state_evolution" begin
 	N = 2
 	M = 2
 	angles_1 = [0.5, 0] * π
 	angles_2 = [0, 0.25, 0] * π
    angles = [angles_1, angles_2]

    ϕ = 0.0
	wf_coeffs = [cis(2 * n * ϕ * π) for n in 0:N - 1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)

    Ψ_out_numerical = mesh_evolution(Ψ_init, angles)
    Ψ_out_analytical = explicit_state_evolution(Ψ_init, angles)
    @test Ψ_out_numerical ≈ Ψ_out_analytical

    N = 4
	M = 4

	angles_1 = [0.25, 0.25, 0.25, 0.25] * π
	angles_2 = [0.25, 0.25, 0.25, 0.25, 0.25] * π
	angles_3 = [0, 0, 0, 0, 0, 0] * π
	angles_4 = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25] * π
    angles = [angles_1, angles_2, angles_3, angles_4]

    #ϕ = 0.3
	#wf_coeffs = [cis(2 * n * ϕ * π) for n in 0:N - 1]
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

@testset "explicit_fs_projection_sp" begin
	angles_1 = [0.5, 0] * π
	angles_2 = [0, 0.25, 0] * π

	# angles_1 = [0.25, 0.25] * π
	# angles_2 = [0.25, 0.25, 0.25] * π

	angles = [angles_1, angles_2]
    @test explicit_fs_projection_sp(0, 0, angles) == ([], [])
    @test explicit_fs_projection_sp(0, 1, angles) == ([], [])
    @test all(explicit_fs_projection_sp(1, 0, angles) .≈
        ([1, 3], 1 / √2 .* Complex.([-1, 1])))
    @test explicit_fs_projection_sp(1, 1, angles) == ([], [])
    @test explicit_fs_projection_sp(2, 0, angles) == ([], [])
    @test all(explicit_fs_projection_sp(2, 1, angles) .≈
        ([1, 3], -im / √2 .* Complex.([1, 1])))
    @test explicit_fs_projection_sp(3, 0, angles) == ([], [])
end

@testset "explicit_fs_projection" begin
	N = 2
    M = 2

    # angles_1 = [0.5, 0] * π
	# angles_2 = [0, 0.25, 0] * π

	angles_1 = [0.25, 0.25] * π
	angles_2 = [0.25, 0.25, 0.25] * π
	angles = [angles_1, angles_2]

    j = lcmk2j(N + M, 0, 0, 0, 0)
    @test all(explicit_fs_projection(j, angles) .≈
        ([1], 1 / 4 .* Complex.([1])))

    j = lcmk2j(N + M, 1, 0, 1, 0)
    @test all(explicit_fs_projection(j, angles) .≈
        ([1, 3, 9, 11], 1 / 4 .* Complex.([1, -1, -1, 1])))

    j = lcmk2j(N + M, 2, 1, 2, 1)
    @test all(explicit_fs_projection(j, angles) .≈
        ([1, 3, 9, 11], 1 / 4 .* Complex.([-1, -1, -1, -1])))

    j = lcmk2j(N + M, 2, 1, 0, 0)
    @test all(explicit_fs_projection(j, angles) .≈
        ([1, 9], -im / 4 .* Complex.([1, 1])))

    N = 4
    M = 2 * (N - 1)

    angles = angles_single_setup(N)
    j_idx_arr_contr_all_symb = Vector{Int64}[]
    j_idx_arr_contr_all_mesh = Vector{Int64}[]
    coeff_arr_all_symb = Vector{ComplexF64}[]
    coeff_arr_all_mesh = Vector{ComplexF64}[]

    for j in 1:N_LOOPS2 * (N + M)^2
        j_idx_arr_contr_symb, coeff_arr_symb =
            TimeBinEncoding._explicit_fs_projection_symbolic_backend(N, M, j, angles)
        push!(j_idx_arr_contr_all_symb, j_idx_arr_contr_symb)
        push!(coeff_arr_all_symb, coeff_arr_symb)

        j_idx_arr_contr_mesh, coeff_arr_mesh =
            TimeBinEncoding._explicit_fs_projection_mesh_backend(N, M, j, angles)
        push!(j_idx_arr_contr_all_mesh, j_idx_arr_contr_mesh)
        push!(coeff_arr_all_mesh, coeff_arr_mesh)
    #=     if !(coeff_arr_mesh ≈ coeff_arr_symb)
            println("j: ", j)
            println("j_idx_arr_contr_symb: ", j_idx_arr_contr_symb)
            println("j_idx_arr_contr_mesh: ", j_idx_arr_contr_mesh)
            println("coeff_arr_symb: ", coeff_arr_symb)
            println("coeff_arr_mesh: ", coeff_arr_mesh)
        end =#
    end

    @test all(j_idx_arr_contr_all_symb .== j_idx_arr_contr_all_mesh)
    #println("coeff_arr_all_symb: ", coeff_arr_all_symb)
    #println("coeff_arr_all_mesh: ", coeff_arr_all_mesh)

   # println("coeff_arr_all_mesh: ", coeff_arr_all_mesh - coeff_arr_all_symb)

    @test all(coeff_arr_all_symb .≈ coeff_arr_all_mesh)
end
