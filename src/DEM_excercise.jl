module DEM_excercise
using StaticArrays
using LinearAlgebra
using TestItems

# Write your package code here.
include("Particles.jl")
include("System.jl")
include("sorting.jl")
export DEMParticles, ParticleCellList
end
