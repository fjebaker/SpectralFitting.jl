using Test
using SpectralFitting

model = PhotoelectricAbsorption() * PowerLaw(a = FitParam(4.0, frozen = true))
cache = make_parameter_cache(model)
SpectralFitting.update_free_parameters!(cache, [2.0, 0.5])
@test cache.parameters == [2.0, 4.0, 0.5]
SpectralFitting.update_frozen_parameters!(cache, [50.0])
@test cache.parameters == [2.0, 50.0, 0.5]
