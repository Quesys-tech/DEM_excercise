struct DEMSimulation{D,T} 
    p::DEMParticles{T}
    max_𝐱::SVector{D,T}
    min_𝐱::SVector{D,T}
    
    step::Int
    Δt::T
    t_end::T
end