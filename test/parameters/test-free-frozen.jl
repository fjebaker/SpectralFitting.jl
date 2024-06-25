using Test
using SpectralFitting

include("../dummies.jl")

model = DummyMultiplicative() * PowerLaw(a = FitParam(4.0, frozen = true))
cache = make_parameter_cache(model)
SpectralFitting.update_free_parameters!(cache, [2.0, 0.5])
@test cache.parameters == [2.0, 4.0, 0.5, 5.0]


# TODO: test to check that auto diff gradients with various sizes can be propagated correctly
