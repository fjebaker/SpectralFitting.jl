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

@test parameter_count(model) == 3
@test free_parameter_count(model) == 2
@test frozen_parameter_count(model) == 1

model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())

@test parameter_count(model) == 8
@test free_parameter_count(model) == 5
@test frozen_parameter_count(model) == 3

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

# test a slightly more complex composite model
cm = DummyMultiplicative() * (DummyAdditive() + DummyAdditiveWithManyFrozen())

params = freeparameters(cm)
expected = [FitParam(1.0), FitParam(2.0), FitParam(1.0), FitParam(1.0), FitParam(1.0)]
@test isapprox.(params, expected) |> all

params = frozenparameters(cm)
expected = [
    FitParam(1.0),
    FitParam(5.0),
    FitParam(2.0),
    FitParam(2.0),
    FitParam(2.0),
    FitParam(2.0),
    FitParam(2.0),
    FitParam(2.0),
    FitParam(2.0),
    FitParam(5.0),
    FitParam(5.0),
]


# test updating parameters on single model
model = DummyAdditive()
new_model = SpectralFitting.updateparameters(model, K = FitParam(2.0))
@test modelparameters(new_model) == [FitParam(2.0), FitParam(1.0), FitParam(5.0)]

new_model = SpectralFitting.updateparameters(model, K = FitParam(2.0), b = FitParam(1.0))
@test modelparameters(new_model) == [FitParam(2.0), FitParam(1.0), FitParam(1.0)]

# test updating parameters on composite model
