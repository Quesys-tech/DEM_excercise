module DEM_excercise
using StaticArrays
using LinearAlgebra
using Tullio
using TestItems

# Write your package code here.
include("Particles.jl")
export DEMParticles
include("CellList.jl")
export ParallelCellList
include("System.jl")
end
