function cell_id(𝐱ᵢ, cl::ParticleCellList{D, T}) where {D, T}
    cell_id = @. floor(Int, (𝐱ᵢ - cl.𝐱_min) / cl.𝐥) + 1
    cartesian = CartesianIndex(cell_id...)
    linear = LinearIndices(cl.max_id)[cartesian]
    return linear
end

function sort!(p::DEMParticles{T}, cell::ParticleCellList{D, T}) where {D, T}
    for i in eachindex(p)
        p[i]._cell_id = cell_id(p[i].𝐱, cell)
    end
    sortperm!(p._permutation, p._cell_id)
    # apply permutation to particles
    p.id = p.id[p._permutation]
    p.type = p.type[p._permutation]
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

function update!(cl::ParticleCellList{D, T}, p::DEMParticles{T}) where {D, T}
    cl.min_id .= typemax(typeof(cl.min_id))
    cl.max_id .= typemin(typeof(cl.max_id))
    for i in eachindex(p)
        cl.min_id[p._cell_id[i]] = min(cl.min_id[p._cell_id[i]], i)
        cl.max_id[p._cell_id[i]] = max(cl.max_id[p._cell_id[i]], i)
    end
end