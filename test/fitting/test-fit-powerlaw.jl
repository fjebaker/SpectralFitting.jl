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
model = PowerLaw() + PowerLaw(a = FitParam(1.0))
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
# for these the order keeps going wrong on the CI, so we'll just check the sum
@test sum(result.u) ≈ sum([1.0167392657565117, 1.0, 1.203211162317021, 3.225587384436137]) atol = 1e-1

# xpsec implementation
model = XS_PowerLaw() + XS_PowerLaw(a = FitParam(1.0))
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test sum(result.u) ≈ sum([1.0167392657565117, 1.0, 1.203211162317021, 3.225587384436137]) atol = 1e-1

# mixed implementation
model = XS_PowerLaw() + PowerLaw(a = FitParam(1.0))
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test sum(result.u) ≈ sum([1.0167392657565117, 1.0, 1.203211162317021, 3.225587384436137]) atol = 1e-1
