
@testset "beam_splitter_operator" begin
    @test beam_splitter_operator(0) ≈ [[1,0] [0,1]]
    @test beam_splitter_operator(π/2) ≈ [[0,im] [im,0]]
    @test beam_splitter_operator(π/4) ≈ 1/√2 *[[1,im] [im,1]]
end