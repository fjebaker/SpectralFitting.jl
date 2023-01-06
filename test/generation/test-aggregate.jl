using Test
using SpectralFitting

include("../dummies.jl")
model = DummyAdditive()

info = SpectralFitting.FunctionGeneration.assemble_aggregate_info(typeof(model), Float64)
@test length(info.infos) == 1
@test length(info.closure_params) == 0
@test info.models == Type[typeof(model)]
@test info.maximum_flux_count == 1
@test length(info.statements) == 1

# composite additive models
new_model = model + model

info =
    SpectralFitting.FunctionGeneration.assemble_aggregate_info(typeof(new_model), Float64)
@test length(info.infos) == 2
@test length(info.closure_params) == 0
@test info.models == Type[typeof(model), typeof(model)]
@test info.maximum_flux_count == 2
@test length(info.statements) == 3


# errors when using two mutliplicatives, would seem to get the parameter parsing wrong
model = DummyMultiplicative() * DummyMultiplicative() * (DummyAdditive() + DummyAdditive())

info = SpectralFitting.FunctionGeneration.assemble_aggregate_info(typeof(model), Float64)