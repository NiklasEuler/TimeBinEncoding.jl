module TimeBinEncoding

using ArgCheck
using LaTeXStrings
using LinearAlgebra
using Random
using SparseArrays
using StatsBase

# Write your package code here.
include("BeamSplitterAngles.jl") # Generation & modification of beam splitter configurations
include("CoherenceExtraction.jl") # Extraction of coherence information from measured data
include("Constants.jl") # Numerical constants
include("ExplicitEvolution.jl") # Explicit evolution based on analytical computations
include("IndexConversions.jl") # Basis index conversion tools
include("MeshEvolution.jl") # Numerical evolution of quantum states in the mesh lattice
include("StateUtility.jl") # Setup and tools for quantum wave functions and density matrices
include("SymbolicEvolution.jl") # Symbolic computing tools for time evolution
include("Visualization.jl") # Visualization for explicit and symbolic results

include("identical_photons/IndexConversionsIdentical.jl")
# Basis index conversion tools for identical photons
include("identical_photons/VisualizationIdentical.jl")
# Visualization for identical photons
include("identical_photons/CoherenceExtractionIdentical.jl")
# Extraction of coherence information from measured data for identical photons
include("identical_photons/ExplicitEvolutionIdentical.jl")
# Explicit evolution based on analytical computations for identical photons
include("identical_photons/MeshEvolutionIdentical.jl")
# Numerical evolution of quantum states in the mesh lattice for identical photons
include("identical_photons/StateUtilityIdentical.jl")
# Setup and tools for quantum wave functions and density matrices for identical photons


end
