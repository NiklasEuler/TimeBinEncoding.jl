
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

#=  @testset "angles_phase_estimation" begin
    N = 4
    angles_12_1 = [0.5, 0, 0, 0] * π
    angles_12_2 = [0, 0.25, 0, 0, 0] * π
    angles_12 = [angles_12_1, angles_12_2]

    angles_23_1 = [0, 0.5, 0, 0] * π
    angles_23_2 = [0, 0, 0.25, 0, 0] * π
    angles_23 = [angles_23_1, angles_23_2]

    angles_34_1 = [0, 0, 0.5, 0] * π
    angles_34_2 = [0, 0, 0, 0.25, 0] * π
    angles_34 = [angles_34_1, angles_34_2]

    angles_all = [angles_12, angles_23, angles_34]

    @test all(angles_all .≈ angles_phase_estimation(N))
 end =#

@testset "angles_compound" begin
    N = 4
    angles_compound_all = angles_compound(N)
#    angles_k1 = angles_phase_estimation(N)

    angles_14_1 = [0.5, 0, 0, 0] * π
    angles_14_2 = [0, 0, 0, 0, 0] * π
    angles_14_3 = [0, 0.5, 0, 0, 0, 0] * π
    angles_14_4 = [0, 0, 0.25, 0.25, 0, 0, 0] * π
    angles_k3 = [angles_14_1, angles_14_2, angles_14_3, angles_14_4]

    @test all(angles_compound_all[1] .≈ angles_k3)

    angles_13_1 = [0.5, 0.5, 0, 0] * π
    angles_13_2 = [0, 0, 0, 0, 0] * π
    angles_13_3 = [0, 0, 0.25, 0.25, 0, 0] * π

    angles_13 = [angles_13_1, angles_13_2, angles_13_3]

    angles_k2 = angles_13

    @test all(angles_compound_all[2] .≈ angles_k2)

    angles_12_1 = [0.5, 0, 0.5, 0] * π
    angles_12_2 = [0, 0.25, 0, 0.25, 0] * π

    angles_k1 = [angles_12_1, angles_12_2]

    @test all(angles_compound_all[3] .≈ angles_k1)
end

#= @testset "angles_compound" begin
    N = 4
    angles_compound_all = angles_compound(N)
#    angles_k1 = angles_phase_estimation(N)

    angles_12_1 = [0.5, 0, 0, 0] * π
    angles_12_2 = [0, 0.25, 0, 0, 0] * π

    angles_12 = [angles_12_1, angles_12_2]

    angles_23_1 = [0, 0.5, 0, 0] * π
    angles_23_2 = [0, 0, 0.25, 0, 0] * π

    angles_23 = [angles_23_1, angles_23_2]

    angles_34_1 = [0, 0, 0.5, 0] * π
    angles_34_2 = [0, 0, 0, 0.25, 0] * π

    angles_34 = [angles_34_1, angles_34_2]

    angles_k1 = [angles_12, angles_23, angles_34]

    @test all(angles_compound_all[1] .≈ angles_k1)

    angles_13_1 = [0.5, 0, 0, 0] * π
    angles_13_2 = [0, 0, 0, 0, 0] * π
    angles_13_3 = [0, 0, 0.25, 0, 0, 0] * π

    angles_13 = [angles_13_1, angles_13_2, angles_13_3]

    angles_24_1 = [0, 0.5, 0, 0] * π
    angles_24_2 = [0, 0, 0, 0, 0] * π
    angles_24_3 = [0, 0, 0, 0.25, 0, 0] * π

    angles_24 = [angles_24_1, angles_24_2, angles_24_3]

    angles_k2 = [angles_13, angles_24]

    @test all(angles_compound_all[2] .≈ angles_k2)

    angles_14_1 = [0.5, 0, 0, 0] * π
    angles_14_2 = [0, 0, 0, 0, 0] * π
    angles_14_3 = [0, 0, 0, 0, 0, 0] * π
    angles_14_4 = [0, 0, 0, 0.25, 0, 0, 0] * π
    angles_14 = [angles_14_1, angles_14_2, angles_14_3, angles_14_4]

    angles_k3 = [angles_14]

    @test all(angles_compound_all[3] .≈ angles_k3)
end =#

@testset "angles_phase_estimation" begin
    N = 4
    ϵ = 0.0
    angles = angles_phase_estimation(N)

    @test all(angles .≈ [
            [0.25, 0.25, 0.25, 0.25],
            [0, 0.25, 0.25, 0.25, 0],
        ] .* π
    )
    @test all(angles_phase_estimation(N, ϵ) .≈ angles)
end


@testset "angles4bins" begin
    N = 8
    l = 0
    m = 3
    p = 5
    q = 6

    angles4bins_all = angles4bins(N, l, m, p, q)

    angles_4b_1 = [
        [0.5, 0, 0, 0, 0, 0.5, 0, 0],
        [0, 0, 0, 0, 0, 0, 0.25, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0],
    ] * π
    @test all(angles4bins_all[1] .≈ angles_4b_1)

    angles_4b_2 = [
        [0.5, 0, 0, 0.5, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0],
    ] * π
    @test all(angles4bins_all[2] .≈ angles_4b_2)

    angles_4b_3 = [
        [0.5, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0.25, 0.25, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0],
    ] * π
    @test (all(angles4bins_all[3] ≈ angles_4b_3))

end

@testset "angles_ref_bin_all_pairs" begin
    N = 2
    angles_ref = angles_ref_bin_all_pairs(N, 0)
    @test angles_ref == [[0.5, 0], [0, 0.25, 0]] * π
    @test angles_ref == angles_ref_bin_all_pairs(N, 1)
    @test_throws ArgumentError angles_ref_bin_all_pairs(N, 2)

    angles_40 = [
        [Θ_13, 0, 0, 0],
        [0.25, 0, 0, 0, 0],
        [0.5, 0, 0, 0, 0, 0],
        [0, 0.25, 0.25, 0.25, 0, 0, 0],
        ] * π
    angles_ref = angles_ref_bin_all_pairs(4, 0)
    @test all(angles_ref .≈ angles_40)

    angles_40_pop = [
        [Θ_13, θ_pop_ref, θ_pop_ref, θ_pop_ref],
        [0.25, 0, 0, 0, 0],
        [0.5, 0, 0, 0, 0, 0],
        [0, 0.25, 0.25, 0.25, 0, 0, 0],
        ] * π
    angles_ref = angles_ref_bin_all_pairs(4, 0; population_bins=true)
    @test all(angles_ref .≈ angles_40_pop)


    angles_41 = [
        [0, Θ_13, 0, 0],
        [0.5, 0.25, 0, 0, 0],
        [0, 0.25, 0.25, 0.25, 0, 0],
        ] * π
    angles_ref = angles_ref_bin_all_pairs(4, 1)
    @test all(angles_ref .≈ angles_41)


    angles_42 = [
        [0.5, 0.5, Θ_23, 0.5],
        [0, 0, 0, 0.25, 0.5],
        [0, 0, 0.25, 0.25, 0.25, 0],
        ] * π
    angles_ref = angles_ref_bin_all_pairs(4, 2)
    @test all(angles_ref .≈ angles_42)

    angles_42_pop = [
        [θ_pop_ref, θ_pop_ref, Θ_23, 0.5],
        [0, 0, 0, 0.25, θ_pop_ref],
        [0, 0, 0.25, 0.25, 0.25, 0],
        ] * π
    angles_ref = angles_ref_bin_all_pairs(4, 2, population_bins=true)
    @test all(angles_ref .≈ angles_42_pop)

    angles_43 = [
        [0.5, 0.5, 0.5, Θ_23],
        [0, 0, 0, 0, 0.25],
        [0, 0, 0, 0, 0, 0.5],
        [0, 0, 0, 0.25, 0.25, 0.25, 0],
        ] * π
    angles_ref = angles_ref_bin_all_pairs(4, 3)
    @test all(angles_ref .≈ angles_43)

    angles_84 = [
        [0, 0, 0, 0, asin(sqrt(4 / 7)) / π, 0.5, 0.5, 0.5],
        [0, 0, 0, 0, Θ_13, asin(sqrt(1 / 4)) / π, 0, 0, 0],
        [0.5, 0.5, 0.5, 0.5, 0.25, 0, Θ_13, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0.25, 0.5, 0.5, 0.5],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0, 0, 0],
        ] * π
        angles_ref = angles_ref_bin_all_pairs(8, 4)
        @test all(angles_ref .≈ angles_84)
end

@testset "graph_coloring" begin
    N = 4
    pairings = graph_coloring(N)

    @test length(pairings) == N - 1
    @test length(pairings[1]) == N / 2
    @test length(unique(reduce(hcat, pairings))) == N * (N - 1) / 2

    perfect_matches_sort = [sort(reduce(vcat, perf_match)) for perf_match in pairings]
    @test all([match == 0:N-1 for match in perfect_matches_sort])

    N = 7
    pairings_odd = graph_coloring(N)

    @test length(pairings_odd) == N
    @test length(pairings_odd[1]) == (N - 1) / 2
    @test length(unique(reduce(hcat, pairings_odd))) == N * (N - 1) / 2

    reduced_pairings = reduce(hcat,(reduce(hcat, pairings_odd)))
    @test N - 1 ∈ reduced_pairings
    @test N ∉ reduced_pairings

    perfect_matches_length =
        [length(unique(reduce(vcat, perf_match))) for perf_match in pairings_odd]
    @test all(perfect_matches_length .== N - 1)
end

@testset "angles_pairs_from_coloring" begin
    N = 4
    pairings = graph_coloring(N)

    angles_1 = angles_pairs_from_coloring(N, pairings[1])
    @test all(
        angles_1 .≈ [
            [0.5, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0.5, 0, 0, 0, 0],
            [0, 0, 0.25, 0.25, 0, 0, 0],
        ] * π
    )

    angles_2 = angles_pairs_from_coloring(N, pairings[2])
    @test all(
        angles_2 .≈ [
            [0.5, 0.5, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0.25, 0.25, 0, 0],
        ] * π
    )

    angles_3 = angles_pairs_from_coloring(N, pairings[3])
    @test all(
        angles_3 .≈ [
            [0.5, 0, 0.5, 0],
            [0, 0.25, 0, 0.25, 0],
        ] * π
    )

    N = 6
    pairings = graph_coloring(N)

    angles_1 = angles_pairs_from_coloring(N, pairings[1]; population_bins=true)
    @test all(
        angles_1 .≈ [
            [θ_pop_ref, 0, 0, θ_pop_ref, θ_pop_ref, θ_pop_ref],
            [0, 0, 0, 0, 0, 0, 0],
            [0, θ_pop_ref, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, θ_pop_ref, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0.25, 0.25, 0.25, 0, 0, 0, 0, 0],
        ] * π
    )

    angles_2 = angles_pairs_from_coloring(N, pairings[2]; population_bins=true)
    @test all(
        angles_2 .≈ [
            [0, θ_pop_ref, θ_pop_ref, 0, θ_pop_ref, θ_pop_ref],
            [0, 0, 0, 0, 0, 0, 0],
            [θ_pop_ref, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, θ_pop_ref, 0, 0, 0, 0, 0],
            [0, 0, 0.25, 0, 0.25, 0.25, 0, 0, 0, 0],
        ] * π
    )


    angles_5 = angles_pairs_from_coloring(N, pairings[5]; population_bins=true)
    @test all(
        angles_5 .≈ [
            [θ_pop_ref, 0, 0, θ_pop_ref, 0, θ_pop_ref],
            [0, 0, θ_pop_ref, 0, 0, 0, 0],
            [0, θ_pop_ref, 0, 0, θ_pop_ref, 0, 0, 0],
            [0, 0, 0.25, 0.25, 0, 0.25, 0, 0, 0],
        ] * π
    )

end

@testset "angles_kth_neighbor_interference" begin

    N = 4
    k = 2
    i = 2
    angles_k = angles_kth_neighbor_interference(N, k)
    angles_k_i = TimeBinEncoding._angles_kth_neighbor(N, k, i)

    @test all(angles_k[i] .≈ angles_k_i)

    angles_k_i_manual = [[0, 0.5, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0.25, 0, 0]] .* π
    @test all(angles_k_i_manual .≈ angles_k_i)

    ϵ = 0
    angles_k_no_noise = angles_kth_neighbor_interference(N, k, ϵ)
    angles_k_i_no_noise = TimeBinEncoding._angles_kth_neighbor(N, k, i, ϵ)

    @test all(angles_k_no_noise[i] .≈ angles_k_i_no_noise)
    @test all(angles_k_i_manual .≈ angles_k_i_no_noise)

    Random.seed!(4441)
    ϵ = 0.01
    angles_k_noise = angles_kth_neighbor_interference(N, k, ϵ)
    angles_k_i_noise = TimeBinEncoding._angles_kth_neighbor(N, k, i, ϵ)

    @test all(isapprox.(angles_k_i_noise, angles_k_noise[i], atol=5e-2))
    @test all(isapprox.(angles_k_i_manual, angles_k_i_noise, atol=5e-2))
end
