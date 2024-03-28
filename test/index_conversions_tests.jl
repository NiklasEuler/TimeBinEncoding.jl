
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

@testset "lcmk <-> j identical" begin

    N = 7
    d_hilbert_space = N_LOOPS2 * N^2 - N

    j_reconstr = zeros(Int64,d_hilbert_space)
	for i in 1:d_hilbert_space
		idxs = j2lcmk_identical(N, i)
		j_reconstr[i] = lcmk2j_identical(N, idxs...)
	end

    @test 1:(d_hilbert_space) == j_reconstr
    @test j2lcmk_identical(N, j_reconstr[end]) == (N - 1, 1, N - 1, 1)
    @test j2lcmk_identical(N, 1) == (0, 0, 0, 0)
    @test j2lcmk_identical(N, 124) == (4, 1, 0, 1)

    N = 12
    d_hilbert_space = N_LOOPS2 * N^2 - N

    j_reconstr = zeros(Int64, d_hilbert_space)
	for i in 1:d_hilbert_space
		idxs = j2lcmk_identical(N, i)
		j_reconstr[i] = lcmk2j_identical(N, idxs...)
	end

    @test 1:(d_hilbert_space) == j_reconstr
    @test j2lcmk_identical(N, j_reconstr[end]) == (N - 1, 1, N - 1, 1)
    @test j2lcmk_identical(N, 1) == (0, 0, 0, 0)
    @test j2lcmk_identical(N, 170) == (3, 1, 2, 0)

    @test_throws ArgumentError lcmk2j_identical(N, 5, 1, 5, 0)
end
