using Test
using SpectralFitting

include("../dummies.jl")
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
