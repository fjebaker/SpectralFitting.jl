using Test
using SpectralFitting

include("../dummies.jl")
include("../fuzz.jl")

model = DummyAdditive()

params = SpectralFitting.all_parameter_symbols(model)
@test params == (:K, :a, :b)
free_symbols = SpectralFitting.free_parameter_symbols(model)
@test free_symbols == (:K, :a)
frozen_symbols = SpectralFitting.frozen_parameter_symbols(model)
@test frozen_symbols == (:b,)

# model info
info = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
@test info.symbols == [:K, :a, :b]
@test info.free == [:K, :a]
@test info.frozen == [:b]

# model introspetion
T = SpectralFitting.FunctionGeneration.model_T(typeof(model))
@test T === typeof(FitParam(1.0))

T = SpectralFitting.numbertype(model)
@test T === typeof(1.0)

# fuzz with all models we have
for model in FUZZ_ALL_MODELS
    _ = SpectralFitting.all_parameter_symbols(model)
    _ = SpectralFitting.free_parameter_symbols(model)
    _ = SpectralFitting.frozen_parameter_symbols(model)
    _ = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
    TT = SpectralFitting.FunctionGeneration.model_T(typeof(model))
    @test TT === typeof(FitParam(1.0))
    TT = SpectralFitting.numbertype(model)
    @test TT === typeof(1.0)
end


# introspection for composite models
cm = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())

@test_throws "" SpectralFitting.all_parameter_symbols(cm)
@test_throws "" SpectralFitting.free_parameter_symbols(cm)
@test_throws "" SpectralFitting.frozen_parameter_symbols(cm)

T = SpectralFitting.FunctionGeneration.model_T(typeof(cm))
@test T == typeof(FitParam(1.0))
T = SpectralFitting.numbertype(cm)
@test T == typeof(1.0)

info = SpectralFitting.FunctionGeneration.getinfo(typeof(cm))
@test length(info) == 3

# order of parameters when we query for them is important
model = DummyAdditiveWithManyFrozen()
res = SpectralFitting.FunctionGeneration.all_parameter_symbols(typeof(model))
@test res == (:K, :a, :b, :c, :d, :e, :f, :g, :h, :i, :j)

res = SpectralFitting.FunctionGeneration.free_parameter_symbols(typeof(model))
@test res == (:K, :h)

res = SpectralFitting.FunctionGeneration.frozen_parameter_symbols(typeof(model))
@test res == (:a, :b, :c, :d, :e, :f, :g, :i, :j)
