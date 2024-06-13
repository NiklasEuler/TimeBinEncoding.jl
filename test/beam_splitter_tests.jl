
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

 @testset "angles_phase_estimation" begin
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
 end

@testset "angles_compound" begin
    N = 4
    angles_compound_all = angles_compound(N)
    angles_k1 = angles_phase_estimation(N)
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
end

@testset "angles_phase_estimation" begin
    N = 4
    ϵ = 0.0
    @test all(angles_phase_estimation(N, ϵ) .≈ angles_phase_estimation(N))
end
