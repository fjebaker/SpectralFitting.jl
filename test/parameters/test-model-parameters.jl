using Test
using SpectralFitting

include("../dummies.jl")
include("../fuzz.jl")

model = DummyAdditive()

params_tuple = @inferred SpectralFitting.parameter_tuple(model)
expected = (FitParam(1.0), FitParam(1.0), FitParam(5.0))
@test isapprox.(params_tuple, expected) |> all

model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())

# test but now composite model
params = SpectralFitting.parameter_tuple(model)
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
@test isapprox.(params, expected) |> all
