@testset "combined_projector_weights_auto" begin
    N = 4
    projector_weights = combined_projector_weights_auto(N)
    @test projector_weights == [2.0, -2.0, 4.0, 1.0, 0.0]
    projector_weights = combined_projector_weights_auto(N, [0, 0, 1, 0, 0])
    @test projector_weights == [0, 0, 0, 4 ,0]
end

@testset "j_out_hom" begin
    N = 4
    j_out_hom_arr = j_out_hom(N, 1, 2)
    @test j_out_hom_arr == [556, 1111, 667]
    @test_throws ArgumentError j_out_hom(N, 3, 3)
    @test_throws ArgumentError j_out_hom(N, 2, 1)
    @test_throws ArgumentError j_out_hom(N, 2, 4)
end

@testset "coherence_extraction_identical" begin
    N = 4
    ϵ = 0.1
    ϵ_angles = 0.05 * π

    Random.seed!(4441)

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

    contr_j_tuples=correlated_short_bins_tuples_identical(N; extract_diagonal=false)

    l = 0
	m = 2
	angles_distance_k = angles_kth_neighbor_interference(N, m - l)
	angles_lm = angles_distance_k[l + 1]
	N_post_lm = N + length(angles_lm)
    j_out_lm_hom = [lcmk2j_super_identical(N_post_lm, m, 0, m+1, 1, m, 0, m+1, 1)]

    j1 = lcmk2j_super_identical(N, l, 0, l, 0, l, 0, l, 0)
    j2 = lcmk2j_super_identical(N, m, 0, m, 0, m, 0, m, 0)
    j_unphys = lcmk2j_super_identical(N, l, 0, l, 0, m, 0, m, 0)
    j_fake = lcmk2j_super_identical(N, m+1, 1, m+1, 1, m+1, 1, m+1, 1)


    projector_weights = [1]

    contr_j_tuples_fake = [(j1, j_fake), (j_fake, j1)]
    contr_j_tuples = [(j1, j2), (j2, j1)]
    contr_j_tuples_some_missing = [(j1, j2), (j2, j1), (j1, j_fake), (j_fake, j1)]
    contr_j_tuples_unphys = [(j1, j2), (j2, j1), (j1, j_unphys), (j_unphys, j1)]

    pop_fs = explicit_fs_pop_identical(ρ_pure, j_out_lm_hom, angles_lm, projector_weights)

    @test_throws ArgumentError coherence_extraction_identical(
        N,
        j_out_lm_hom,
        pops_pure,
        pop_fs,
        angles_lm,
        contr_j_tuples_fake,
        projector_weights
    )

    coh_extr = coherence_extraction_identical(
        N,
        j_out_lm_hom,
        pops_pure,
        pop_fs,
        angles_lm,
        contr_j_tuples,
        projector_weights
    )

    @test coh_extr ≈ ρ_pure[j2, j1] + ρ_pure[j1, j2]

    @test_logs (Logging.Warn,"Some of the scheduled coherences have a vanishing weight in the "*
        "given final-state projectors. Please check again and consider adapting the sche"*
        "duled coherences in `contr_j_tuples`."
    ) coherence_extraction_identical(
        N,
        j_out_lm_hom,
        pops_pure,
        pop_fs,
        angles_lm,
        contr_j_tuples_some_missing,
        projector_weights
    )

    @test_logs (Logging.Warn,"The coherences scheduled for extraction have differing weights in "*
    "the chosen final-state projectors and can thus not straight-forwardly be extrac"*
    "ted. Please check input again."
    ) coherence_extraction_identical(
        N,
        j_out_lm_hom,
        pops_pure,
        pop_fs,
        angles_lm,
        contr_j_tuples_unphys,
        projector_weights
    )

end

@testset "combined_measurement_coherence_extraction_identical" begin
    N = 4
    ϵ = 0.1
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

    projector_weights = [1, 0, 0, 0, 0]

    combined_weights_4b, pop_fs_4b =
        combined_weights_pops_4bins_all(N, ρ_pure, projector_weights, phases=true)
	combined_weights_hom, pop_fs_hom =
        combined_weights_pops_hom(N, ρ_pure, projector_weights)

    combined_weights = combined_weights_hom + combined_weights_4b
	pop_fs_combined = pop_fs_4b + pop_fs_hom

    @test_logs (Logging.Warn,"The coherences scheduled for extraction have differing weights in "*
        "the chosen final-state projectors and can thus not straight-forwardly be extrac"*
        "ted. Please check input again."
        ) coh_extract_reconstr = combined_measurement_coherence_extraction_identical(
            N,
            combined_weights,
            pops_pure,
            pop_fs_combined,
            contr_j_tuples
        )

    projector_weights = [0, 0, 1, 0, 0]

    combined_weights_4b, pop_fs_4b =
        combined_weights_pops_4bins_all(N, ρ_pure, projector_weights, phases=true)
	combined_weights_hom, pop_fs_hom =
        combined_weights_pops_hom(N, ρ_pure, projector_weights)

    combined_weights = combined_weights_hom + combined_weights_4b
	pop_fs_combined = pop_fs_4b + pop_fs_hom

    @test_logs (Logging.Warn,"Some of the scheduled coherences have a vanishing weight in the "*
        "given final-state projectors. Please check again and consider adapting the sche"*
        "duled coherences in `contr_j_tuples`."
    ) (Logging.Warn,"The coherences scheduled for extraction have differing weights in "*
    "the chosen final-state projectors and can thus not straight-forwardly be extrac"*
    "ted. Please check input again."
    ) coh_extract_reconstr = combined_measurement_coherence_extraction_identical(
            N,
            combined_weights,
            pops_pure,
            pop_fs_combined,
            contr_j_tuples
    )
end
