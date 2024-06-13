@testset "lcmk <-> j identical" begin

    N = 7
    d_hilbert_space = Int(N_LOOPS * N * (N_LOOPS * N + 1) / 2)

    j_reconstr = zeros(Int64,d_hilbert_space)
	for i in 1:d_hilbert_space
		idxs = j2lcmk_identical(N, i)
		j_reconstr[i] = lcmk2j_identical(N, idxs...)
	end

    @test 1:(d_hilbert_space) == j_reconstr
    @test j2lcmk_identical(N, j_reconstr[end]) == (N - 1, 1, N - 1, 1)
    @test j2lcmk_identical(N, 1) == (0, 0, 0, 0)
    @test j2lcmk_identical(N, 61) == (2, 1, 2, 1)

    N = 12
    d_hilbert_space = Int(N_LOOPS * N * (N_LOOPS * N + 1) / 2)

    j_reconstr = zeros(Int64, d_hilbert_space)
	for i in 1:d_hilbert_space
		idxs = j2lcmk_identical(N, i)
		j_reconstr[i] = lcmk2j_identical(N, idxs...)
	end

    @test 1:(d_hilbert_space) == j_reconstr
    @test j2lcmk_identical(N, j_reconstr[end]) == (N - 1, 1, N - 1, 1)
    @test j2lcmk_identical(N, 1) == (0, 0, 0, 0)
    @test j2lcmk_identical(N, 70) == (1, 1, 1, 1)

    @test_throws ArgumentError lcmk2j_identical(N, 5, 1, 5, 0)
end
