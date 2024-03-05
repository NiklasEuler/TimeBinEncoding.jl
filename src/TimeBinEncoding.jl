module TimeBinEncoding

using ArgCheck
using LaTeXStrings
using LinearAlgebra
using Random
using SparseArrays
using StatsBase

# Write your package code here.
include("MeshEvolution.jl")
include("SymbolicEvolution.jl")
include("ExplicitEvolution.jl")
include("IndexConversions.jl")
include("Visualization.jl")
include("StatePrep.jl")
include("CoherenceExtraction.jl")

end
