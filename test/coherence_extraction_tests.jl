
@testset "coherence_extraction_tests" begin
	N = 4
	M = 4
    ϵ = 0.1
    ϵ_angles = 0.05 * π

    ϕ = 0
    wf_coeffs = [cis.(2 * ϕ * k * π) for k in 1:N]
	Ψ_init = insert_initial_state(correlated_timebin_state(wf_coeffs))
	ρ_pure = density_matrix(Ψ_init)
    pops_pure = populations(ρ_pure)

    Ψ_mes =  insert_initial_state(correlated_timebin_state(fill(1 / sqrt(N), N)))
	mes_fidelity = fidelity(Ψ_mes, ρ_pure)

	j_out_arr = [lcmk2j(N + M, 3, 0, 3, 0), lcmk2j(N + M, 4, 1, 4, 1)]
	angles_1 = [0.5, 0.5, 0, 0] * π
	angles_2 = [0, 0, 0, 0, 0] * π
	angles_3 = [0, 0, 0.25, 0.25, 0, 0] * π
	angles_4 = [0, 0, 0.25, 0.25, 0.25, 0, 0] * π
	angles = [angles_1, angles_2, angles_3, angles_4]

    #coherence_extraction(N, j_out_arr, ρ_pure, angles)

    pop_fs = explicit_fs_projection_expval(ρ_pure, j_out_arr, angles)

    @test isapprox((@test_logs min_level=Logging.Warn coherence_extraction(N, j_out_arr,
    pops_pure, pop_fs, angles)), mes_fidelity, atol = 1e-8)

	angles_1 = [0.2, 0.6, 0.1, 0] * π
	angles_2 = [0, 0.8, 0.2, 0, 0] * π
	angles_3 = [0, 0.35, 0.41, 0.9, 0, 0] * π
	angles_4 = [0, 0, 0.12, 0.26, 0.83, 0, 0] * π
	angles = [angles_1, angles_2, angles_3, angles_4]
    pop_fs = explicit_fs_projection_expval(ρ_pure, j_out_arr, angles)
	@test_throws ArgumentError coherence_extraction(N, j_out_arr, pops_pure, pop_fs, angles)

	N = 4
	M = 2

    ϕ = 0
	wf_coeffs = [cis(n * ϕ * π) for n in 0:N - 1]
	tb_state = correlated_timebin_state(wf_coeffs)
	Ψ_init = insert_initial_state(tb_state)
	ρ_pure = density_matrix(Ψ_init)
    pops_pure = populations(ρ_pure)

	angles_1_1 = [0.5, 0, 0.5, 0] * π
	angles_1_2 = [0, 0.25, 0, 0.25, 0] * π
	angles_1 = [angles_1_1, angles_1_2]
	j_01 = [lcmk2j(N + M, 1, 0, 1, 0), lcmk2j(N + M, 2, 1, 2, 1)]

	extract_diagonal = false
    pop_fs = explicit_fs_projection_expval(ρ_pure, j_01, angles_1)

	@test isapprox((@test_logs (:warn, "Some of the scheduled coherences have a vanishing "*
        "weight in the given final-state projectors. Please check again and consider "*
        "adapting the scheduled coherences in `contr_j_idxs`.") min_level =
        Logging.Warn coherence_extraction(N, j_01, pops_pure, pop_fs, angles_1;
        extract_diagonal=extract_diagonal)), 1 / 8, atol = 1e-8)
end

@testset "compound coherence extraction" begin
    N = 8
    ϵ = 0.1
    ϵ_angles = 0.05 * π
    ϕ = 0.00

    wf_coeffs = [cis(n * ϕ * π) for n in 0:N - 1]
    Ψ_init = insert_initial_state(correlated_timebin_state(wf_coeffs))
	ρ_pure = density_matrix(Ψ_init)
    pops_pure = populations(ρ_pure)

    Ψ_mes =  insert_initial_state(correlated_timebin_state(fill(1 / sqrt(N), N)))
	mes_fidelity = fidelity(Ψ_mes, ρ_pure)

    angles_compound_all = angles_compound(N)
    pops_fs_all_pure = pops_fs_compound(ρ_pure, angles_compound_all)
    @test isapprox(
        (@test_logs min_level=Logging.Warn compound_coherence_extraction(
            pops_pure, pops_fs_all_pure
        )),
        mes_fidelity,
        atol = 1e-8
    )
    ρ_mixed = density_matrix_dephased(Ψ_init, ϵ)
    pops_mixed = populations(ρ_mixed)

    pops_fs_all_mixed = pops_fs_compound(ρ_mixed, angles_compound_all)

    angles_compound_all_noisy = angles_compound(N, ϵ_angles)
    pops_fs_all_pure_noisy = pops_fs_compound(ρ_pure, angles_compound_all_noisy)

    @test compound_coherence_extraction(pops_mixed, pops_fs_all_mixed) ≤
        compound_coherence_extraction(pops_pure, pops_fs_all_pure)
    @test compound_coherence_extraction(pops_pure, pops_fs_all_pure_noisy) ≤
        compound_coherence_extraction(pops_pure, pops_fs_all_pure)
end

@testset "angles_single_setup" begin
    N = 4
    angles_1_mod = [0.5, 0.5, 0, 0] * π
	angles_2_mod = [0, 0, 0, 0, 0] * π
	angles_3_mod = [0, 0, 0.25, 0.25, 0, 0] * π
	angles_4_mod = [0, 0, 0.5, 0.25, 0.5, 0, 0] * π
	angles_5_mod = [0, 0, 0, 0, 0, 0, 0, 0] * π
	angles_6_mod = [0, 0, 0, 0, 0.25, 0, 0, 0, 0] * π
	angles_mod = [angles_1_mod, angles_2_mod, angles_3_mod, angles_4_mod, angles_5_mod,
        angles_6_mod]
    angles = angles_single_setup(N)
    @test angles == angles_mod

    N = 8
    M_mod = 14
	angles_8_mod = [zeros(Float64, n) for n in N:N + M_mod - 1]
	angles_8_mod[1][1:4] .= 0.5 * π
	angles_8_mod[5][5:8] .= 0.25 * π
	angles_8_mod[7][[5, 10]] .= 0.5 * π
	angles_8_mod[7][7:8] .= 0.25 * π
	angles_8_mod[8][[6, 8, 10]] .= 0.25 * π
	angles_8_mod[8][[7, 9]] .= 0.5 * π
	angles_8_mod[9][[6, 11]] .= 0.5 * π
	angles_8_mod[10][[9]] .= 0.25 * π
	angles_8_mod[12][10] = 0.25 * π
	angles_8_mod[14][11] = 0.25 * π
    angles_8 = angles_single_setup(N)
    @test angles_8 == angles_8_mod

    N = 16
    M_mod = 30
    angles_16_mod_symm = [zeros(Float64, n) for n in N:N + M_mod - 1]
	angles_16_mod_symm[1][1:8] .= 0.5 * π
	angles_16_mod_symm[9][9:16] .= 0.25 * π
	angles_16_mod_symm[13][[9, 10, 19, 20]] .= 0.5 * π
	angles_16_mod_symm[13][13:16] .= 0.25 * π
	angles_16_mod_symm[15][[13, 18]] .= 0.5 * π
	angles_16_mod_symm[15][[11, 12, 15, 16, 19, 20]] .= 0.25 * π
	angles_16_mod_symm[16][[11, 13, 15, 17, 19, 21]] .= 0.5 * π
	angles_16_mod_symm[16][[12, 14, 16, 18, 20]] .= 0.25 * π
	angles_16_mod_symm[17][[14, 19]] .= 0.5 * π
	angles_16_mod_symm[18][[13, 17, 21]] .= 0.25 * π
	angles_16_mod_symm[19][[12, 13, 22, 23]] .= 0.5 * π
	angles_16_mod_symm[20][[18]] .= 0.25 * π
	angles_16_mod_symm[22][[19]] .= 0.25 * π
	angles_16_mod_symm[24][[20]] .= 0.25 * π
	angles_16_mod_symm[26][[21]] .= 0.25 * π
	angles_16_mod_symm[28][[22]] .= 0.25 * π
	angles_16_mod_symm[30][[23]] .= 0.25 * π
    angles_16 = angles_single_setup(N)
    @test angles_16 == angles_16_mod_symm

    N = 6
    @test_throws ArgumentError angles_single_setup(N)
end

@testset "j_out_single_setup" begin
    N = 8
    M_mod = 14
    j_mod_arr = [lcmk2j(N + M_mod, 7, 0, 7, 0), lcmk2j(N + M_mod, 14, 1, 14, 1),
        lcmk2j(N + M_mod, 8, 0, 8, 0), lcmk2j(N + M_mod, 13, 1, 13, 1),
        lcmk2j(N + M_mod, 9, 0, 9, 0), lcmk2j(N + M_mod, 12, 1, 12, 1),
        lcmk2j(N + M_mod, 10, 0, 10, 0), lcmk2j(N + M_mod, 11, 1, 11, 1),
    ]
    j_arr = j_out_single_setup(N)
    @test all([j in j_arr for j in j_mod_arr])
    @test all([j in j_mod_arr for j in j_arr])

    N = 10
    @test_throws ArgumentError j_out_single_setup(N)
end
