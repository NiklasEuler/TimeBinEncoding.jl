module TimeBinEncoding

using ArgCheck
using LaTeXStrings
using LinearAlgebra
using Random
using SparseArrays
using StatsBase

# Write your package code here.
include("BeamSplitterAngles.jl") # Generation & modification of beam splitter configurations
include("Constants.jl") # Numerical constants
include("MeshEvolution.jl") # Numerical evolution of quantum states in the mesh lattice
include("SymbolicEvolution.jl") # Symbolic computing tools for time evolution
include("ExplicitEvolution.jl") # Explicit evolution based on analytical computations
include("IndexConversions.jl") # Basis index conversion tools
include("Visualization.jl") # Visualization for explicit and symbolic results
include("StatePrep.jl") # Setup and readout for quantum wave functions and density matrices
include("CoherenceExtraction.jl") # Extraction of coherence information from measured data

end
