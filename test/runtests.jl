using Test
using SpectralFitting

include("utils.jl")

@testset "function-generation" verbose = true begin
    @testset "aggregation" begin
        include("generation/test-aggregate.jl")
    end
    @testset "utilities" begin
        include("generation/test-parsing-utilities.jl")
    end
end

@testset "macro" verbose = true begin
    @testset "xspecmodel" begin
        include("macros/test-xspec.jl")
    end
end

@testset "composite-algebra" verbose = true begin
    @testset "model-algebra" begin
        include("composite/test-algebra.jl")
    end
    @testset "model-invocation" begin
        include("composite/test-invocation.jl")
    end
end

@testset "parameters" verbose = true begin
    @testset "model-parameters" begin
        include("parameters/test-model-parameters.jl")
    end
    @testset "fit-params" begin
        include("parameters/test-fit-params.jl")
    end
end

@testset "model-library" verbose = true begin
    @testset "table-models" begin
        include("models/test-table-models.jl")
    end
    @testset "julia-models" begin
        include("models/test-julia-models.jl")
    end

    # only test XSPEC models when not using CI
    # since model data access is annoying
    @ciskip @testset "xspec-models" begin
        include("models/test-xspec-models.jl")
    end
    @testset "general-models" begin
        include("models/test-general-models.jl")
        # include the general xspec models only when not CI
        @ciskip include("models/test-general-xspec-models.jl")
    end

    @testset "model-consistency" begin
        include("models/test-model-consistency.jl")
    end
end

@testset "io" verbose = true begin
    @testset "printing" begin
        include("io/test-printing.jl")
    end
end

@testset "fitting" verbose = true begin
    @testset "powerlaws" begin
        include("fitting/test-fit-powerlaw.jl")
    end 
end

using Aqua

# little bit of aqua
Aqua.test_undefined_exports(SpectralFitting)
Aqua.test_unbound_args(SpectralFitting)
