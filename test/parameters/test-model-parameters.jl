using Test
using SpectralFitting

include("../dummies.jl")
include("../fuzz.jl")

model = DummyAdditive()

# test but single model
params = modelparameters(model)
expected = [FitParam(1.0), FitParam(1.0), FitParam(5.0)]
@test isapprox.(params, expected) |> all

params_tuple = @inferred SpectralFitting.model_parameters_tuple(model)
expected = (FitParam(1.0), FitParam(1.0), FitParam(5.0))
@test isapprox.(params_tuple, expected) |> all

@test parameter_count(model) == 3

model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())

@test parameter_count(model) == 8

# test but now composite model
params = modelparameters(model)
expected = [
    FitParam(1.0),
    FitParam(1.0),
    FitParam(5.0),
    FitParam(1.0),
    FitParam(1.0),
    FitParam(5.0),
    FitParam(1.0),
    FitParam(5.0),
]
@test isapprox.(params, expected) |> all

params_tuple = SpectralFitting.model_parameters_tuple(model)
expected = (
    FitParam(1.0),
    FitParam(1.0),
    FitParam(5.0),
    FitParam(1.0),
    FitParam(1.0),
    FitParam(5.0),
    FitParam(1.0),
    FitParam(5.0),
)
@test isapprox.(params_tuple, expected) |> all

# test updating parameters on single model
model = DummyAdditive()
new_model = SpectralFitting.updateparameters(model, K = FitParam(2.0))
@test modelparameters(new_model) == [FitParam(2.0), FitParam(1.0), FitParam(5.0)]

new_model = SpectralFitting.updateparameters(model, K = FitParam(2.0), b = FitParam(1.0))
@test modelparameters(new_model) == [FitParam(2.0), FitParam(1.0), FitParam(1.0)]
