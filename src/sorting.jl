function cell_id(𝐱ᵢ, cl::CellList{D, T}) where {D, T}
    cell_id = @. floor(Int, (𝐱ᵢ - cl.bc.𝐱_min) / cl.h) + 1
    cartesian = CartesianIndex(cell_id...)
    linear = LinearIndices(cl.id_max)[cartesian]
    return linear
end

function sort!(p::DEMParticles{T}, cell::CellList{BC, D, T}) where {BC, D, T}
    Threads.@threads for i in 1:length(p)
        p._cell_id[i] = cell_id(p.𝐱[:, i], cell)
    end
    sortperm!(p._permutation, p._cell_id)
    # apply permutation to particles
    p.id .= p.id[p._permutation]
    p.type .= p.type[p._permutation]
    p.𝐱 .= p.𝐱[:, p._permutation]
    p.𝐮 .= p.𝐮[:, p._permutation]
    p.𝛚 .= p.𝛚[:, p._permutation]
    p.𝛉 .= p.𝛉[:, p._permutation]
    p.m .= p.m[p._permutation]
    p.Π .= p.Π[p._permutation]
    p.r .= p.r[p._permutation]
    p._cell_id .= p._cell_id[p._permutation]
    p._sorted = true
    return p
end
@testitem "sort!" begin
    using LinearAlgebra
    import DEM_excercise: DEMParticles, CellList, sort!

    𝐱 = [7.4 6.9 7.4 0.8 1.0 5.4 9.5 0.1 4.4 3.5
         3.8 9.9 2.1 1.4 0.4 6.6 2.7 1.4 6.2 0.3
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    N = size(𝐱, 2)
    type = ones(Int, N)
    𝛉 = zero(𝐱)
    𝐮 = zero(𝐱)
    𝛚 = zero(𝐱)
    p = DEMParticles(type, rand(N), rand(N), rand(N) .* 2, 𝐱, 𝛉, 𝐮, 𝛚)
    bc = NoneBoundaryCondition([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
    cl = CellList(2.0, bc)
    sort!(p, cl)
    @test p.𝐱[:, 1] == [0.8, 1.4, 0]
    @test p.𝐱[:, 10] == [6.9, 9.9, 0]
end

function update!(cl::CellList{BC, D, T}, p::DEMParticles{T}) where {BC, D, T}
    sort!(p, cl)
    cl.id_min .= -1
    cl.id_max .= -1

    cl.id_min[p._cell_id[1]] = 1
    Threads.@threads for i in (1 + 1):(length(p) - 1)
        ibl = p._cell_id[i - 1]
        ibi = p._cell_id[i]
        ibr = p._cell_id[i + 1]

        ibl < ibi && (cl.id_min[ibi] = i)
        ibr > ibi && (cl.id_max[ibi] = i)
    end
    cl.id_min[p._cell_id[length(p)]] = length(p)
end
@testitem "update!" begin
    using LinearAlgebra
    import DEM_excercise: DEMParticles, CellList, sort!, update!

    𝐱 = [7.4 6.9 7.4 0.8 1.0 5.4 9.5 0.1 4.4 3.5
         3.8 9.9 2.1 1.4 0.4 6.6 2.7 1.4 6.2 0.3
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    N = size(𝐱, 2)
    type = ones(Int, N)
    𝛉 = zero(𝐱)
    𝐮 = zero(𝐱)
    𝛚 = zero(𝐱)
    p = DEMParticles(type, rand(N), rand(N), rand(N) .* 2, 𝐱, 𝛉, 𝐮, 𝛚)
    bc = NoneBoundaryCondition([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
    cl = CellList(2.0, bc)
    update!(cl, p)
    @test cl.id_min[1, 1, 1] == 1
    @test cl.id_max[1, 1, 1] == 3
    @test cl.id_min[2, 1, 1] == 4
    @test cl.id_max[2, 1, 1] == 4
    @test cl.id_min[3, 1, 1] == -1
    @test cl.id_min[3, 1, 1] == -1
end
