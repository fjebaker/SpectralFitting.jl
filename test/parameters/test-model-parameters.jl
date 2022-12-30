using Test
using SpectralFitting

include("../dummies.jl")
include("../fuzz.jl")

model = DummyAdditive()

# test but single model
params = modelparameters(model)
expected = [FitParam(1.0), FitParam(1.0), FitParam(5.0)]
@test isapprox.(params, expected) |> all

params = freeparameters(model)
expected = [FitParam(1.0), FitParam(1.0)]
@test isapprox.(params, expected) |> all

params = frozenparameters(model)
expected = [FitParam(5.0)]
@test isapprox.(params, expected) |> all


model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())

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

params = freeparameters(model)
expected = [FitParam(1.0), FitParam(1.0), FitParam(1.0), FitParam(1.0), FitParam(1.0)]
@test isapprox.(params, expected) |> all

params = frozenparameters(model)
expected = [FitParam(5.0), FitParam(5.0), FitParam(5.0)]
@test isapprox.(params, expected) |> all
