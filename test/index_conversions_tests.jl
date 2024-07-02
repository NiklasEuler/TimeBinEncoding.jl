
@testset "lcmk <-> j" begin
    N = 5
    j_arr = 1:(4*N^2)
    lcmk_arr = zeros(4, 4*N^2)
    j_arr_reconstructed = zeros(4*N^2)
    for j in eachindex(j_arr)
       lcmk = collect(j2lcmk(N, j))
       lcmk_arr[:, j] = lcmk
       j_arr_reconstructed[j] = lcmk2j(N, lcmk...)
    end
    @test j_arr == j_arr_reconstructed
    @test lcmk_arr[:, 1] == [0, 0, 0, 0]
    @test lcmk_arr[:, end] == [N - 1, 1, N - 1, 1]
    @test lcmk_arr[:, 77] == [3, 1, 3, 0]
    @test_throws ArgumentError j2lcmk(0, 1)
    @test_throws ArgumentError j2lcmk(N, 0)
    @test_throws ArgumentError j2lcmk(N, 4*N^2 + 1)
    @test_throws ArgumentError lcmk2j(N, 4, 0, N, 0)
    @test_throws ArgumentError lcmk2j(N, N, 0, 2, 0)
    @test_throws ArgumentError lcmk2j(N, 4, 0, -1, 0)
    @test_throws ArgumentError lcmk2j(N,-2, 0, 2, 0)
    @test_throws ArgumentError lcmk2j(N, 4, 100, 2, 0)
    @test_throws ArgumentError lcmk2j(N, 1, 0, 2, 500)
    @test_throws ArgumentError lcmk2j(N, 3, 0, 2, -1)
    @test_throws ArgumentError lcmk2j(N, 2, -1, 2, 1)
end

@testset "lm <-> j" begin
    N = 5
    j_arr = 1:(N^2)
    lm_arr = zeros(2, N^2)
    j_arr_reconstructed = zeros(N^2)
    for j in eachindex(j_arr)
       lm = collect(j2lm(N, j))
       lm_arr[:, j] = lm
       j_arr_reconstructed[j] = lm2j(N, lm...)
    end
    @test j_arr == j_arr_reconstructed
    @test lm_arr[:, 1] == [0, 0]
    @test lm_arr[:, end] == [N - 1, N - 1]
    @test lm_arr[:, 13] == [2, 2]
    @test_throws ArgumentError j2lm(0, 1)
    @test_throws ArgumentError j2lm(N, 0)
    @test_throws ArgumentError j2lm(N, N^2 + 1)
    @test_throws ArgumentError lm2j(N, 4, N)
    @test_throws ArgumentError lm2j(N, N, 2)
    @test_throws ArgumentError lm2j(N, 4, -1)
    @test_throws ArgumentError lm2j(N,-2, 2)
end

@testset "lc <-> j" begin
    N = 5
    j_arr = 1:(2 * N)
    lc_arr = zeros(2, 2 * N)
    j_arr_reconstructed = zeros(2 * N)
    for j in eachindex(j_arr)
       lc = collect(j2lc(j))
       lc_arr[:, j] = lc
       j_arr_reconstructed[j] = lc2j(lc...)
    end
    @test j_arr == j_arr_reconstructed
    @test lc_arr[:, 1] == [0, 0]
    @test lc_arr[:, end] == [N - 1, 1]
    @test lc_arr[:, 6] == [2, 1]
    @test_throws ArgumentError j2lc(0)
    @test_throws ArgumentError lc2j(-1, 1)
    @test_throws ArgumentError lc2j(1, -1)
    @test_throws ArgumentError lc2j(4, 2)
end


@testset "correlated_short_bins_tuples" begin
    N = 4
    tuples = correlated_short_bins_tuples(N)
    @test length(tuples) == N * (N - 1)
    j1 = lcmk2j(N, 0, 0, 0, 0)
    j2 = lcmk2j(N, 1, 0, 1, 0)
    @test (j1, j2) in tuples
    @test (j2, j1) in tuples
    @test (j1, j1) ∉ tuples

    N = 5

    tuples = correlated_short_bins_tuples(N; extract_diagonal = true)
    @test length(tuples) == N^2

    j2 = lcmk2j(N, N - 1, 0, N - 1, 0)
    @test (j2, j2) in tuples
end
