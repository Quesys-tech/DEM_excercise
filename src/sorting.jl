function cell_id(𝐱ᵢ, cl::CellList{D, T}) where {D, T}
    cell_id = @. floor(Int, (𝐱ᵢ - cl.𝐱_min) / cl.𝐥) + 1
    cartesian = CartesianIndex(cell_id...)
    linear = LinearIndices(cl.max_id)[cartesian]
    return linear
end

function sort!(p::DEMParticles{T}, cell::CellList{D, T}) where {D, T}
    for i in eachindex(p)
        p[i]._cell_id = cell_id(p[i].𝐱, cell)
    end
    sortperm!(p._permutation, p._cell_id)
    # apply permutation to particles
    p.id .= p.id[p._permutation]
    p.type .= p.type[p._permutation]
    p.𝐱 .= p.𝐱[:, p._permutation]
    p.𝐮 .= p.𝐮[:, p._permutation]
    p.𝛚 .= p.𝛚[:, p._permutation]
    p.𝛉 .= p.𝛉[:, p._permutation]
    p.m = p.m[p._permutation]
    p.Π = p.Π[p._permutation]
    p.r = p.r[p._permutation]
    p._cell_id = p._cell_id[p._permutation]
    p._sorted = true

    return p
end
@testitem "sort!" begin
    using LinearAlgebra
    import DEM_excercise: DEMParticles, CellList, sort!
    
    𝐱 = [7.4 6.9 7.4 0.8 1.0 5.4 9.5 0.1 4.4 3.5
         3.8 9.9 2.1 1.4 0.4 6.6 2.7 1.4 6.2 0.3
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    p = DEMParticles(type, rand(10), rand(10), rand(10).*2, 𝐱, 𝛉, 𝐮, 𝛚)
    sort!(p, cl)
end

function update!(cl::CellList{D, T}, p::DEMParticles{T}) where {D, T}
    cl.min_id .= typemax(typeof(cl.min_id))
    cl.max_id .= typemin(typeof(cl.max_id))
    for i in eachindex(p)
        cl.min_id[p._cell_id[i]] = min(cl.min_id[p._cell_id[i]], i)
        cl.max_id[p._cell_id[i]] = max(cl.max_id[p._cell_id[i]], i)
    end
end