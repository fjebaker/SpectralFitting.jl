using Test
using SpectralFitting

@testset "composite-algebra" begin
    include("composite/test-algebra.jl")
end

@testset "generation" begin
    include("generation/test-aggregate.jl")
    include("generation/test-parsing-utilities.jl")
end