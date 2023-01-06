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

params = freeparameters(model)
expected = [FitParam(1.0), FitParam(1.0)]
@test isapprox.(params, expected) |> all

params_tuple = @inferred SpectralFitting.free_parameters_tuple(model)
expected = (FitParam(1.0), FitParam(1.0))
@test isapprox.(params_tuple, expected) |> all

params = frozenparameters(model)
expected = [FitParam(5.0)]
@test isapprox.(params, expected) |> all

params_tuple = @inferred SpectralFitting.frozen_parameters_tuple(model)
expected = (FitParam(5.0),)
@test isapprox.(params_tuple, expected) |> all

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
model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())
# new_model = SpectralFitting.updateparameters(model, K_1 = FitParam(2.0))
# modelparameters(new_model)


# TODO: model rebuilding
# SpectralFitting.FunctionGeneration.rebuild_composite_model(typeof(model))
# SpectralFitting.FunctionGeneration.getinfo(typeof(model))
