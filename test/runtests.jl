using DEM_excercise
using Test
using Aqua

@testset "DEM_excercise.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(DEM_excercise)
    end
    # Write your tests here.
end
