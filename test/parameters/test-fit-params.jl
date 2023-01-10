using Test
using SpectralFitting

include("../utils.jl")
include("../dummies.jl")

# test that fit parameters can be updated in models easily
model = DummyAdditive()
updateparameters!(model; K = 2.0)
@test model.K.value == 2.0

updateparameters!(model; K = 4.0, b = 3.0)
@test model.K.value == 4.0
@test model.b.value == 3.0

updateparameters!(model; a = FitParam(0.0, lower_limit = -100.0))
@test model.a.value == 0.0
@test model.a.lower_limit == -100.0

# test that same interface works for composite models
