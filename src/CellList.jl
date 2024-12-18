abstract type AbstractBoundaryCondition{D, T <: Real} end
struct NoneBoundaryCondition{D, T} <: AbstractBoundaryCondition{D, T}
    𝐱_min::SVector{D, T}
    𝐱_max::SVector{D, T}
    function NoneBoundaryCondition(
            𝐱_min::AbstractVector{T}, 𝐱_max::AbstractVector{T}) where {T}
        @assert length(𝐱_min) == length(𝐱_max)
        D = length(𝐱_min)
        new{D, T}(SVector{D}(𝐱_min), SVector{D}(𝐱_max))
    end
end
@testitem "NoneBoundaryCondition" begin
    using StaticArrays
    bc = NoneBoundaryCondition([0.0, 0.0], [1.0, 1.0])
    @test bc.𝐱_min == SVector(0.0, 0.0)
    @test bc.𝐱_max == SVector(1.0, 1.0)
end
struct PeriodicBoundaryCondition{D, T} <: AbstractBoundaryCondition{D, T}
    𝐱_min::SVector{D, T}
    𝐱_max::SVector{D, T}
    function PeriodicBoundaryCondition(
            𝐱_min::AbstractVector{T}, 𝐱_max::AbstractVector{T}) where {T}
        @assert length(𝐱_min) == length(𝐱_max)
        D = length(𝐱_min)
        new{D, T}(SVector{D}(𝐱_min), SVector{D}(𝐱_max))
    end
end
@testitem "PeriodicBoundaryCondition" begin
    using StaticArrays
    bc = PeriodicBoundaryCondition([0.0, 0.0], [1.0, 1.0])
    @test bc.𝐱_min == SVector(0.0, 0.0)
    @test bc.𝐱_max == SVector(1.0, 1.0)
end

abstract type AbstractCellList{BC <: AbstractBoundaryCondition, D, T} end
struct CellList{BC, D, T} <: AbstractCellList{BC, D, T}
    h::T
    bc::BC
    id_min::Array{Int, D}
    id_max::Array{Int, D}
    function CellList(
            h::T, bc::BC) where {BC <: AbstractBoundaryCondition{D, T}} where {D, T}
        num_cell_dim = @. ceil(Int, (bc.𝐱_max - bc.𝐱_min) / h)
        id_min = zeros(Int, num_cell_dim...)
        id_max = zeros(Int, num_cell_dim...)
        new{BC, D, T}(h, bc, id_min, id_max)
    end
end

@testitem "CellList" begin
    using StaticArrays
    bc = PeriodicBoundaryCondition([0.0, 0.0], [1.0, 1.0])
    cl = CellList(0.1, bc)
    @test cl.h == 0.1
    @test cl.id_min == zeros(Int, 10, 10)
    @test cl.id_max == zeros(Int, 10, 10)
end
