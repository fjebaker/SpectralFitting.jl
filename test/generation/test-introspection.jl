using Test
using SpectralFitting

include("../dummies.jl")


model = DummyAdditive()
infos = model_parameter_info(model)
@test infos == (
    (:K, FitParam(1.0), true),
    (:a, FitParam(1.0), true),
    (:b, FitParam(5.0), false)
)

model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())
infos = model_parameter_info(model)
@test infos == [
    (:K, FitParam(1.0), true),
    (:a, FitParam(1.0), true),
    (:b, FitParam(5.0), false),
    (:K, FitParam(1.0), true),
    (:a, FitParam(1.0), true),
    (:b, FitParam(5.0), false),
    (:a, FitParam(1.0), true),
    (:b, FitParam(5.0), false),
]

model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())
