using Test
using SpectralFitting

testdir = get(
    ENV,
    "SF_TEST_SUITE_DATA",
    @__DIR__() * "/../../../spectral-fitting-test-suite/sample-data",
)
@show testdir

has_test_dir = isdir(testdir)
if !has_test_dir
    @warn "No test data found. Skipping some tests."
end
if !has_test_dir && get(ENV, "CI", false)
    error("Missing test dir on CI")
end

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
        include("parameters/test-free-frozen.jl")
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
    if has_test_dir
        @testset "datasets" begin
            include("datasets/test-ogip.jl")
            include("datasets/test-grouping.jl")
            include("datasets/test-binning.jl")
            include("datasets/test-datasets.jl")
        end
    else
        @warn "Skipping dataset tests."
    end
end

@testset "fitting" verbose = true begin
    include("fitting/test-fit-simple-dataset.jl")
    include("fitting/test-cache.jl")
    @testset "powerlaws" begin
        include("fitting/test-fit-powerlaw.jl")
    end
    @testset "multifits" begin
        include("fitting/test-fit-multi.jl")
        include("fitting/test-fit-optim.jl")
    end
    if has_test_dir
        @testset "sample-data" begin
            include("fitting/test-sample-data.jl")
        end
    else
        @warn "Skipping dataset tests."
    end
end

using Aqua

# little bit of aqua
Aqua.test_undefined_exports(SpectralFitting)
Aqua.test_unbound_args(SpectralFitting)
