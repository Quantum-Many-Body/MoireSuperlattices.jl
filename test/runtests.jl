using Test
using SafeTestsets

@testset "all" begin
    @time @safetestset "lattices" begin include("lattices.jl") end
    @time @safetestset "systems" begin include("systems.jl") end
    @time @safetestset "analysis" begin include("analysis.jl") end
end
