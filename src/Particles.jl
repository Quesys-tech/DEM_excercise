
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
    m::Vector{T}
    Π::Vector{T}
    r::Vector{T}
    𝐱::Matrix{T}
    𝛉::Matrix{T}
    𝐮::Matrix{T}
    𝛚::Matrix{T}
    _permutation::Vector{Int}
    _cell_id::Vector{Int}
    _sorted::Bool

    function DEMParticles(type::Vector{Int}, m::Vector{T}, Π::Vector{T}, r::Vector{T},
            𝐱::Matrix{T}, 𝛉::Matrix{T}, 𝐮::Matrix{T}, 𝛚::Matrix{T}) where {T}
        N = length(m)
        @assert N == length(Π) == length(r) == size(𝐱, 2) == size(𝛉, 2) == size(𝐮, 2) ==
                size(𝛚, 2)
        @assert size(𝐱, 1) == size(𝛉, 1) == size(𝐮, 1) == size(𝛚, 1) == 3

        id = collect(1:N)
        _permutation = collect(1:N)
        _cell_id = zeros(Int, N)
        _sorted = false
        new{T}(id, type, m, Π, r, 𝐱, 𝛉, 𝐮, 𝛚, _permutation, _cell_id, _sorted)
    end
end
@testitem "DEMParticles" begin
    type = ones(Int, 10)
    m = rand(10)
    Π = rand(10)
    r = rand(10) * 0.1
    𝐱 = [7.4 6.9 7.4 0.8 1.0 5.4 9.5 0.1 4.4 3.5
         3.8 9.9 2.1 1.4 0.4 6.6 2.7 1.4 6.2 0.3
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    𝛉 = zero(𝐱)
    𝐮 = zero(𝐱)
    𝛚 = zero(𝐱)
    p = DEMParticles(type, m, Π, r, 𝐱, 𝛉, 𝐮, 𝛚)
    @test p.id == collect(1:10)
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