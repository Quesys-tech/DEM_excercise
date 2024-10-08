
abstract type AbstractDEMParticle{T <: AbstractFloat} end
"""
    DEMParticles{T <: AbstractFloat}

A struct to store the information of the particles in the DEM simulation.

# Fields
- `id::Vector{Int}`: The id of the particles.
- `type::Vector{Int}`: The type of the particles.
- `𝛉::Matrix{T}`: The rotation matrix of the particles. size: 3xN
- `m::Vector{T}`: The mass of the particles.
- `Π::Vector{T}`: The inertia of the particles.
- `r::Vector{T}`: The radius of the particles.
- `𝐱::Matrix{T}`: The position of the particles. size: 3xN
- `𝐮::Matrix{T}`: The velocity of the particles. size: 3xN
- `𝛚::Matrix{T}`: The angular velocity of the particles. size: 3xN
"""
struct DEMParticles{T <: AbstractFloat} <: AbstractDEMParticle{T}
    id::Vector{Int}
    type::Vector{Int}
    𝛉::Matrix{T}
    m::Vector{T}
    Π::Vector{T}
    r::Vector{T}
    𝐱::Matrix{T}
    𝐮::Matrix{T}
    𝛚::Matrix{T}
    _permutation::Vector{Int}
    _cell_id::Vector{Int}
    _sorted::Bool
end

struct ParticleCellList{D, T}
    h::T
    𝐱_min::SVector{D, T}
    𝐱_max::SVector{D, T}
    min_id::Array{Int, D}
    max_id::Array{Int, D}
    function ParticleCellList(h::T, 𝐱_min::SVector{D, T}, 𝐱_max::SVector{D, T}) where {D, T}
        num_cell_dim = @. ceil(Int, (𝐱_max - 𝐱_min) / h)
        min_id = zeros(Int, num_cell_dim...)
        max_id = zeros(Int, num_cell_dim...)
        new{D, T}(h, 𝐱_min, 𝐱_max, min_id, max_id)
    end
end

@testitem "ParticleCellList" begin
    using StaticArrays
    cl = ParticleCellList(0.1, SVector(0.0, 0.0), SVector(1.0, 1.0))
    @test cl.h == 0.1
    @test cl.𝐱_min == SVector(0.0, 0.0)
    @test cl.𝐱_max == SVector(1.0, 1.0)
    @test cl.min_id == zeros(Int, 10, 10)
    @test cl.max_id == zeros(Int, 10, 10)
end