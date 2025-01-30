@testset "combined_projector_weights_auto" begin
    N = 4
    projector_weights = combined_projector_weights_auto(N)
    @test projector_weights == [2.0, -2.0, 4.0, 1.0, 0.0]
    projector_weights = combined_projector_weights_auto(N, [0, 0, 1, 0, 0])
    @test projector_weights == [0, 0, 0, 4 ,0]
end

@testset "compound coherence extraction" begin
    N = 4
    ϵ = 0.1
    ϵ_angles = 0.05 * π

    Random.seed!(4441)

    d_local_hs_b = Int(N * (N + 1) / 2)
    d_local_hs_bl = N * (2 * N + 1)
    d_full_hs_bl = d_local_hs_bl ^ 2

    Ψ_init = spzeros(ComplexF64, d_full_hs_bl)
	for l in 0:N - 1
		for m in l:N - 1
			j_super = lcmk2j_super_identical(N, l, 0, m, 0, l, 0, m, 0)
            Ψ_init[j_super] = 1
		end
	end
	normalize!(Ψ_init)

	ρ_pure = density_matrix(Ψ_init)
    pops_pure = populations(ρ_pure)

    mes_fidelity = fidelity(Ψ_init, ρ_pure)

    projector_weights = combined_projector_weights_auto(N)

    combined_weights_4b, pop_fs_4b =
        combined_weights_pops_4bins_all(N, ρ_pure, projector_weights, phases=true)
	combined_weights_hom, pop_fs_hom =
        combined_weights_pops_hom(N, ρ_pure, projector_weights)

    combined_weights = combined_weights_hom + combined_weights_4b
	pop_fs_combined = pop_fs_4b + pop_fs_hom

    contr_j_tuples=correlated_short_bins_tuples_identical(N; extract_diagonal=false)

    coh_pop = sum(pops_pure[correlated_short_bins_idxs_identical(N)]) / d_local_hs_b
    coh_extract_true = (mes_fidelity - coh_pop) * d_local_hs_b

	@test isapprox(
        (@test_logs min_level=
            Logging.Warn combined_measurement_coherence_extraction_identical(
                N,
                combined_weights,
                pops_pure,
                pop_fs_combined,
                contr_j_tuples
            )),
            coh_extract_true,
            atol = 1e-8
    )


    ρ_mixed = density_matrix_dephased_krauss_identical(N, ρ_pure, ϵ)
    pops_mixed = populations(ρ_mixed)

    mes_fidelity_mixed = fidelity(Ψ_init, ρ_mixed)

    combined_weights_4b_mixed, pop_fs_4b_mixed =
        combined_weights_pops_4bins_all(N, ρ_mixed, projector_weights, phases=true)
	combined_weights_hom_mixed, pop_fs_hom_mixed =
        combined_weights_pops_hom(N, ρ_mixed, projector_weights)

    combined_weights_mixed = combined_weights_hom_mixed + combined_weights_4b_mixed
	pop_fs_combined_mixed = pop_fs_4b_mixed + pop_fs_hom_mixed

    @test combined_weights_4b_mixed ≤ combined_weights_4b
    @test combined_weights_hom_mixed ≤ combined_weights_hom

    coh_pop_mixed = sum(pops_mixed[correlated_short_bins_idxs_identical(N)]) / d_local_hs_b

    coh_extract_true_mixed = (mes_fidelity_mixed - coh_pop_mixed) * d_local_hs_b

    @test coh_extract_true_mixed ≤ coh_extract_true

    coh_extract_reconstr_mixed = combined_measurement_coherence_extraction_identical(
        N,
        combined_weights_mixed,
        pops_mixed,
        pop_fs_combined_mixed,
        contr_j_tuples
    )

	@test coh_extract_reconstr_mixed ≤ coh_extract_true_mixed
end
