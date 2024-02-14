module TimeBinEncoding

using SparseArrays
using LinearAlgebra
using Random
using ArgCheck

const global n_loops = 2 # number of fiber loops. Saved as const to avoid magic numbers.
const global n_loops2 = 4 # number of fiber loops squared. Saved as const to avoid magic numbers.

# Write your package code here.
include("MeshEvolution.jl")
include("SymbolicEvolution.jl")
include("IndexConversions.jl")

end
