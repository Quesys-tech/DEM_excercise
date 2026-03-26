struct ParallelCellList{D, T}
    domain_size::SVector{D, T}
    cell_id::Vector{Int} # cell id for each particle
    cache_permutation::Vector{Int} # cache for particle permutation
    cutoff::SVector{D, T} # cutoff distance for each dimension
    cell_min_id::Array{Int, D} # min cell id for each cell
    cell_max_id::Array{Int, D} # max cell id for each cell
    neighbor_table::Matrix{Int} # neighbor table for each particle
    num_j_gt_i::Vector{Int} # number of particles j > i for each particle i
    num_j_lt_i::Vector{Int} # number of particles j < i for each particle i
    sum_j_gt_i::Vector{Int} # cumulative sum of num_j_gt_i
    sum_j_lt_i::Vector{Int} # cumulative sum of num_j_lt_i
    pair_i_p::Vector{Int} # particle i in pair list
    pair_j_p::Vector{Int} # particle j in pair list
    reference_table::Matrix{Int} # reference table for each particle
    pair_rel_pos::Matrix{T} # relative position of each pair
    cache_scalar::Vector{T} # cache for scalar values
    cache_vector::Matrix{T} # cache for vector values
    function ParallelCellList(domain_size::SVector{D, T},
            cutoff_in::SVector{D, T},
            num_particles::Integer,
            max_neighbors_in::Integer
    ) where {D, T}
        @assert length(domain_size) == D
        @assert all(domain_size .> cutoff_in .> 0)
        @assert num_particles > 0
        @assert max_neighbors_in > 0
        num_cells = floor.(Int, domain_size ./ cutoff_in)
        cutoff = domain_size ./ num_cells
        max_neighbors = ceil(Int, max_neighbors_in * prod(cutoff) / prod(cutoff_in))
        max_pairs = num_particles * max_neighbors

        cell_id = Vector{Int}(undef, num_particles)
        cache_permutation = Vector{Int}(undef, num_particles)
        cell_min_id = Array{Int, D}(undef, num_cells...)
        cell_max_id = Array{Int, D}(undef, num_cells...)
        neighbor_table = Matrix{Int}(undef, max_neighbors, num_particles)
        num_j_gt_i = Vector{Int}(undef, num_particles)
        num_j_lt_i = Vector{Int}(undef, num_particles)
        sum_j_gt_i = Vector{Int}(undef, num_particles)
        sum_j_lt_i = Vector{Int}(undef, num_particles)
        pair_i_p = Vector{Int}(undef, max_pairs)
        pair_j_p = Vector{Int}(undef, max_pairs)
        reference_table = Matrix{Int}(undef, max_neighbors, num_particles)
        pair_rel_pos = Matrix{T}(undef, D, max_pairs)
        cache_scalar = Vector{T}(undef, max_pairs)
        cache_vector = Matrix{T}(undef, D, max_pairs)
        new{D, T}(domain_size,
            cell_id,
            cache_permutation,
            cutoff,
            cell_min_id,
            cell_max_id,
            neighbor_table,
            num_j_gt_i,
            num_j_lt_i,
            sum_j_gt_i,
            sum_j_lt_i,
            pair_i_p,
            pair_j_p,
            reference_table,
            pair_rel_pos,
            cache_scalar,
            cache_vector)
    end
end

function ParallelCellList(domain_size_in::AbstractVector{T},
        cutoff_in::AbstractVector{T},
        num_particles::Integer,
        max_neighbors_in::Integer
) where {T}
    D = length(domain_size_in)
    @assert length(cutoff_in) == D
    domain_size = SVector{D, T}(domain_size_in...)
    cutoff = SVector{D, T}(cutoff_in...)
    return ParallelCellList(domain_size, cutoff, num_particles, max_neighbors_in)
end

function ParallelCellList(domain_size_in::AbstractVector{T},
        cutoff_in::T,
        num_particles::Integer,
        max_neighbors_in::Integer
) where {T}
    D = length(domain_size_in)
    domain_size = SVector{D, T}(domain_size_in...)
    cutoff = SVector{D, T}(ntuple(_ -> cutoff_in, D)...)
    return ParallelCellList(domain_size, cutoff, num_particles, max_neighbors_in)
end


function cartesian_cell_index(pcl::ParallelCellList{D, T}, position) where {D, T}
    ci = floor.(Int, position ./ pcl.cutoff) .+ 1
    return CartesianIndex(ci...)
end

function cell_index(pcl::ParallelCellList{D, T}, ci::CartesianIndex{D}) where {D, T}
    return LinearIndices(pcl.cell_min_id)[ci]
end
function cell_index(pcl::ParallelCellList{D, T}, position) where {D, T}
    return cell_index(pcl, cartesian_cell_index(pcl, position))
end

function set_cell_id!(
        pcl::ParallelCellList{D, T}, positions::AbstractMatrix{T}) where {D, T}
    @assert size(positions, 2) == length(pcl.cell_id)
    for i in axes(positions, 2)
        @assert all(0 .<= positions[:, i] .< pcl.domain_size)
        pcl.cell_id[i] = cell_index(pcl, positions[:, i])
    end
end

function wrap_position!(
        positions::AbstractMatrix{T}, pcl::ParallelCellList{D, T}) where {D, T}
    @tullio positions[k, i] = mod(positions[k, i], pcl.domain_size[k])
end

function sort_particles!(
        p::DEMParticles{D, T}, pcl::ParallelCellList{D, T}) where {D, T}
    @assert length(p) == length(pcl.cell_id) == length(pcl.cache_permutation)
    wrap_position!(p.𝐱, pcl)
    set_cell_id!(pcl, p.𝐱)
    sortperm!(pcl.cache_permutation, pcl.cell_id)
    permute_particles!(p, pcl.cache_permutation)
    set_cell_id!(pcl, p.𝐱)
end

function update!(pcl::ParallelCellList{D, T}, particles::DEMParticles{D, T}) where {D, T}
    num_particles = length(particles)
    @assert size(pcl.neighbor_table, 2) == num_particles
    @assert length(pcl.num_j_gt_i) == num_particles
    @assert length(pcl.num_j_lt_i) == num_particles
    @assert length(pcl.sum_j_gt_i) == num_particles
    @assert length(pcl.sum_j_lt_i) == num_particles

    sort_particles!(particles, pcl)
    # update cell min/max ids
    @assert size(pcl.cell_min_id) == size(pcl.cell_max_id)
    pcl.cell_min_id .= -1
    pcl.cell_max_id .= -1
    pcl.cell_min_id[pcl.cell_id[1]] = 1
    pcl.cell_max_id[pcl.cell_id[num_particles]] = num_particles
    for i in 2:num_particles
        if pcl.cell_id[i] != pcl.cell_id[i - 1]
            pcl.cell_min_id[pcl.cell_id[i]] = i
        end
        if pcl.cell_id[i - 1] != pcl.cell_id[i]
            pcl.cell_max_id[pcl.cell_id[i - 1]] = i - 1
        end
    end
    # update neighbor table
    pcl.num_j_gt_i .= 0
    pcl.num_j_lt_i .= 0
    pcl.sum_j_gt_i .= 0
    pcl.sum_j_lt_i .= 0
    @assert size(pcl.neighbor_table, 2) == num_particles
    cutoff_inv = one(T) ./ pcl.cutoff

    neighbor_cell_range = CartesianIndices(([-1:1 for k in 1:D]...,))
    for i in 1:num_particles
        pcl.neighbor_table[:, i] .= -1
        i_cell = cartesian_cell_index(pcl, cols[i])
        for offset in neighbor_cell_range
            j_cell_unwrapped = i_cell + offset
            j_cell = CartesianIndex([mod1(j_cell_unwrapped[k], size(pcl.cell_min_id, k))
                                     for k in 1:D]...)
            pcl.cell_max_id[j_cell] < 0 && continue # skip empty cells
            for j in pcl.cell_min_id[j_cell]:pcl.cell_max_id[j_cell]
                i == j && continue
                eta_ij_sq(cols[i], cols[j], pcl.domain_size, cutoff_inv) > one(T) &&
                    continue
                if j > i
                    pcl.num_j_gt_i[i] += 1
                    pcl.neighbor_table[pcl.num_j_gt_i[i], i] = j
                else  # j < i
                    pcl.num_j_lt_i[i] += 1
                    pcl.neighbor_table[end - pcl.num_j_lt_i[i] + 1, i] = j
                end
            end
        end
        if pcl.num_j_gt_i[i] + pcl.num_j_lt_i[i] > size(pcl.neighbor_table, 1)
            throw(ArgumentError("Too many neighbors for particle $i: " *
                                "$(pcl.num_j_gt_i[i] + pcl.num_j_lt_i[i]) > " *
                                "$(size(pcl.neighbor_table, 1))"))
        end
    end
    cumsum!(pcl.sum_j_gt_i, pcl.num_j_gt_i)
    cumsum!(pcl.sum_j_lt_i, pcl.num_j_lt_i)
    # create contact pairs and reference table
    pcl.reference_table .= -1
    for i in 1:num_particles
        for i_table in 1:pcl.num_j_gt_i[i]
            i_list = pcl.sum_j_gt_i[i] - pcl.num_j_gt_i[i] + i_table
            j = pcl.neighbor_table[i_table, i]
            @assert j > i
            @tullio pcl.pair_rel_pos[k, $i_list] = cols[$j][k] - cols[$i][k]
            @tullio pcl.pair_rel_pos[k, $i_list] += -round(
                pcl.pair_rel_pos[k, $i_list] / pcl.domain_size[k]) * pcl.domain_size[k]
            pcl.pair_i_p[i_list] = i
            pcl.pair_j_p[i_list] = j
            pcl.reference_table[i_table, i] = i_list
            for k in 1:pcl.num_j_lt_i[j]
                if pcl.neighbor_table[end - k + 1, j] == i
                    pcl.reference_table[pcl.num_j_gt_i[j] + k, j] = i_list
                end
            end
        end
    end
end

function interaction_symmetric_scalar!(
        func::Function, output::Vector{T}, pcl::ParallelCellList{D, T}) where {D, T}
    for i_list in 1:last(pcl.sum_j_gt_i)
        pcl.cache_scalar[i_list] = func(
            view(pcl.pair_rel_pos, :, i_list), pcl.pair_i_p[i_list], pcl.pair_j_p[i_list])
    end
    for i in eachindex(output)
        output[i] = zero(eltype(output))
        for i_table in 1:pcl.num_j_gt_i[i]
            i_list = pcl.reference_table[i_table, i]
            output[i] += pcl.cache_scalar[i_list]
        end
        for i_table in (pcl.num_j_gt_i[i] + 1):(pcl.num_j_gt_i[i] + pcl.num_j_lt_i[i])
            i_list = pcl.reference_table[i_table, i]
            output[i] += pcl.cache_scalar[i_list]
        end
    end
end

function interaction_antisymmetric_scalar!(
        func::Function, output::Vector{T}, pcl::ParallelCellList{D, T}) where {D, T}
    for i_list in 1:last(pcl.sum_j_gt_i)
        pcl.cache_scalar[i_list] = func(
            view(pcl.pair_rel_pos, :, i_list), pcl.pair_i_p[i_list], pcl.pair_j_p[i_list])
    end
    for i in eachindex(output)
        output[i] = zero(eltype(output))
        for i_table in 1:pcl.num_j_gt_i[i]
            i_list = pcl.reference_table[i_table, i]
            output[i] += pcl.cache_scalar[i_list]
        end
        for i_table in (pcl.num_j_gt_i[i] + 1):(pcl.num_j_gt_i[i] + pcl.num_j_lt_i[i])
            i_list = pcl.reference_table[i_table, i]
            output[i] -= pcl.cache_scalar[i_list]
        end
    end
end

function interaction_antisymmetric_vector!(
        func::Function, output::Matrix{T}, pcl::ParallelCellList{D, T}) where {D, T}
    @assert size(output, 1) == D
    @assert size(output, 2) == size(pcl.reference_table, 2)
    @sync for i_list in 1:last(pcl.sum_j_gt_i)
        Threads.@spawn begin
            view(pcl.cache_vector,
            :,
            i_list) .= func(
                view(pcl.pair_rel_pos, :, i_list), pcl.pair_i_p[i_list], pcl.pair_j_p[i_list])
        end
    end
    @sync for i in axes(output, 2)
        Threads.@spawn begin
            output_i = view(output, :, i)
            output_i .= zero(T)
            for i_table in 1:pcl.num_j_gt_i[i]
                i_list = pcl.reference_table[i_table, i]
                output_i .+= pcl.cache_vector[:, i_list]
            end
            for i_table in (pcl.num_j_gt_i[i] + 1):(pcl.num_j_gt_i[i] + pcl.num_j_lt_i[i])
                i_list = pcl.reference_table[i_table, i]
                output_i .-= pcl.cache_vector[:, i_list]
            end
        end
    end
end
