
@testset "symbolic_fs_projection_sp" begin

    @test_throws ArgumentError symbolic_fs_projection_sp(0, 0, 0)
    @test_throws ArgumentError symbolic_fs_projection_sp(0, 0)

    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs =
        symbolic_fs_projection_sp(1, 0, 0)
    @test j_idx_arr_fs == [1]
    @test trigonometric_history_arr_fs == [fill(0, 1, 1)]
    @test trigonometric_history_arr_fs == [fill(0, 1, 1)]

    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs =
        symbolic_fs_projection_sp(2, 2, 1)
    @test j_idx_arr_fs == [1, 3]
    @test trigonometric_history_arr_fs == [[[1] [0]], [[0] [1]]]
    @test angle_history_arr_fs == [[[0] [1]], [[1] [1]]]

    j = lc2j(2, 1)
    j_idx_arr_fs_j, trigonometric_history_arr_fs_j, angle_history_arr_fs_j =
        symbolic_fs_projection_sp(2, j)
    @test j_idx_arr_fs == j_idx_arr_fs_j
    @test trigonometric_history_arr_fs == trigonometric_history_arr_fs_j
    @test angle_history_arr_fs == angle_history_arr_fs_j

    j_idx_arr_fs, trigonometric_history_arr_fs, angle_history_arr_fs =
        symbolic_fs_projection_sp(2, 0, 1)
    @test j_idx_arr_fs == trigonometric_history_arr_fs == angle_history_arr_fs == []
end
