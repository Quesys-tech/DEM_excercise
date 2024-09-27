
abstract type AbstractDEMParticle{T <: AbstractFloat} end

struct DEMParticles{T <: AbstractFloat} <: AbstractDEMParticle{T}
    id::Vector{Int}
    type::Vector{Int}
    𝐱::Matrix{T}
    𝐮::Matrix{T}
    𝛚::Matrix{T}
    𝛉::Matrix{T}
    m::Vector{T}
    Π::Vector{T}
    r::Vector{T}
    _permutation::Vector{Int}
    _cell_id::Vector{Int}
    _sorted::Bool
end

struct ParticleCellList{D, T <: AbstractFloat}
    𝐥::SVector{D, T}
    𝐱_min::SVector{D, T}
    𝐱_max::SVector{D, T}
    min_id::Array{Int, D}
    max_id::Array{Int, D}
end