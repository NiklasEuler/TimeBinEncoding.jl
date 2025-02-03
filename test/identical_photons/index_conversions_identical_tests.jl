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

@testset "lcmk <-> j super_identical" begin

    N = 7
    d_hilbert_space = Int(N * (2 * N + 1))^2

    j_reconstr = zeros(Int64,d_hilbert_space)
	for i in 1:d_hilbert_space
		idxs = j_super2lcmk_identical(N, i)
		j_reconstr[i] = lcmk2j_super_identical(N, idxs...)
	end

    @test 1:(d_hilbert_space) == j_reconstr
    @test j_super2lcmk_identical(N, j_reconstr[end]) ==
        (N - 1, 1, N - 1, 1, N - 1, 1, N - 1, 1)
    @test j_super2lcmk_identical(N, 1) == (0, 0, 0, 0, 0, 0, 0, 0)
    @test j_super2lcmk_identical(N, 280) == (0, 0, 1, 0, 3, 0, 3, 0)

    N = 12
    d_hilbert_space = Int(N * (2 * N + 1))^2

    j_reconstr = zeros(Int64, d_hilbert_space)
	for i in 1:d_hilbert_space
		idxs = j_super2lcmk_identical(N, i)
		j_reconstr[i] = lcmk2j_super_identical(N, idxs...)
	end

    @test 1:(d_hilbert_space) == j_reconstr
    @test j_super2lcmk_identical(N, j_reconstr[end]) ==
        (N - 1, 1, N - 1, 1, N - 1, 1, N - 1, 1)
    @test j_super2lcmk_identical(N, 1) == (0, 0, 0, 0, 0, 0, 0, 0)
    @test j_super2lcmk_identical(N, 4441) == (0, 0, 7, 0, 6, 1, 9, 1)

    @test_throws ArgumentError lcmk2j_super_identical(N, 5, 1, 5, 0, 1, 0, 2, 0)
    @test_throws ArgumentError lcmk2j_super_identical(N, 2, 1, 3, 0, 4, 0, 2, 0)

    @test_throws ArgumentError j_super2lcmk_identical(N, d_hilbert_space + 1)
end

@testset "correlated_short_bins_tuples" begin

    N = 4

    tuples = correlated_short_bins_tuples_identical(N)
    d_local_hs_b = N * (N + 1) / 2
        # local basis size for two indistinguishable photons, only short loop

    @test length(tuples) == d_local_hs_b * (d_local_hs_b - 1)

    j1 = lcmk2j_super_identical(N, 0, 0, 0, 0, 0, 0, 0, 0)
    j2 = lcmk2j_super_identical(N, 1, 0, 2, 0, 1, 0, 2, 0)
    j3 = lcmk2j_super_identical(N, 3, 0, 3, 0, 3, 0, 3, 0)

    @test (j1, j2) in tuples
    @test (j2, j1) in tuples
    @test (j2, j3) in tuples

    @test (j1, j1) ∉ tuples
    @test (j3, j3) ∉ tuples

    N = 5
    d_local_hs_b = N * (N + 1) / 2
        # local basis size for two indistinguishable photons, only short loop

    tuples = correlated_short_bins_tuples_identical(N; extract_diagonal = true)

    @test length(tuples) == d_local_hs_b^2

    j1 = lcmk2j_super_identical(N, 0, 0, 0, 0, 0, 0, 0, 0)
    j2 = lcmk2j_super_identical(N, N - 1, 0, N - 1, 0, N - 1, 0, N - 1, 0)

    @test (j1, j2) in tuples
    @test (j2, j2) in tuples
    @test (j1, j2) in tuples

end
