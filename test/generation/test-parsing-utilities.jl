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
@test T === typeof(1.0)

T = SpectralFitting.generated_model_parameter_type(model)
@test T === typeof(1.0)

# fuzz with all models we have
for model in FUZZ_ALL_MODELS
    _ = SpectralFitting.all_parameter_symbols(model)
    _ = SpectralFitting.free_parameter_symbols(model)
    _ = SpectralFitting.frozen_parameter_symbols(model)
    _ = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
    TT = SpectralFitting.FunctionGeneration.model_T(typeof(model))
    @test TT === typeof(1.0)
    TT = SpectralFitting.generated_model_parameter_type(model)
    @test TT === typeof(1.0)
end 