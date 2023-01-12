using Test
using SpectralFitting

include("../dummies.jl")

# generate some fake powerlaw data
dummy_data = make_dummy_dataset((E) -> (E^(-3.0)))

# test that both julia and xspec implementations can fit simple
model = PowerLaw()
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test result.u ≈ [1.1824968887784022, 3.0809925004740317] atol = 1e-4
model = XS_PowerLaw()
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test result.u ≈ [1.1824968887784022, 3.0809925004740317] atol = 1e-4

# composite power law
dummy_data = make_dummy_dataset((E) -> (E^(-3.0) + E^(-1)))

# julia implementation
model = PowerLaw() + PowerLaw()
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test result.u ≈
      [1.0929605489222565, 1.0372890824708065, 1.0994740399422531, 3.3105621077461107] atol =
    1e-4

# xpsec implementation
model = XS_PowerLaw() + XS_PowerLaw()
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test result.u ≈
      [1.0929605489222565, 1.0372890824708065, 1.0994740399422531, 3.3105621077461107] atol =
    1e-4

# mixed implementation
model = XS_PowerLaw() + PowerLaw()
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test result.u ≈
      [1.0994740242232883, 3.3105621234900915, 1.092960558240788, 1.0372890866573738] atol =
    1e-4
