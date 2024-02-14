using TimeBinEncoding
using Test

@testset "lcmk <-> j" begin
    # Write your tests here.
    N = 5
    j_arr = 1:(4*N^2)
    lcmk_arr = zeros(4, 4*N^2)
    j_arr_reconstructed = zeros(4*N^2)
    for j in eachindex(j_arr)
       lcmk = collect(j2lcmk(N,j))
       lcmk_arr[:,j] = lcmk
       j_arr_reconstructed[j] = lcmk2j(N,lcmk...)
    end
    @test j_arr == j_arr_reconstructed
    @test lcmk_arr[:,1] == [0,0,0,0]
    @test lcmk_arr[:,end] == [N-1,1,N-1,1]
end
