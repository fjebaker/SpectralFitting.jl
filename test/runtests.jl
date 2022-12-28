using Test
using SpectralFitting

@testset "function-generation" verbose=true begin
    @testset "aggregation" begin
        include("generation/test-aggregate.jl")
    end
    @testset "utilities" begin
        include("generation/test-parsing-utilities.jl")
    end
end

@testset "composite-algebra" verbose=true begin
    @testset "model-algebra" begin
        include("composite/test-algebra.jl")
    end
    @testset "model-invocation" begin
        include("composite/test-invocation.jl")
    end
end

@testset "model-library" verbose=true begin
    @testset "table-models" begin
        include("models/test-table-models.jl")
    end
    @testset "julia-models" begin
        include("models/test-julia-models.jl")
    end
    @testset "xspec-models" begin
        include("models/test-xspec-models.jl")
    end
    @testset "general-models" begin
        include("models/test-models.jl")
    end
end
