using LinearAlgebra
using Logging
using SparseArrays
using Suppressor
using Random
using Test
using TimeBinEncoding

Random.seed!(8675309) # / Jenny

@testset "coherence_extraction_tests" begin
    include("coherence_extraction_tests.jl")
end

@testset "explicit_tests" begin
    include("explicit_tests.jl")
end

@testset "index_conversions_tests" begin
    include("index_conversions_tests.jl")
end

@testset "mesh_tests" begin
    include("mesh_tests.jl")
end

@testset "state_prep_tests" begin
    include("state_prep_tests.jl")
end

@testset "symbolic_tests" begin
    include("symbolic_tests.jl")
end

@testset "visualization_tests" begin
    include("visualization_tests.jl")
end
