module DEM_excercise
using StaticArrays
using LinearAlgebra
using TestItems

# Write your package code here.
include("Particles.jl")
export DEMParticles
include("CellList.jl")
export PeriodicBoundaryCondition, NoneBoundaryCondition, CellList
include("System.jl")
include("sorting.jl")
export update!
end
