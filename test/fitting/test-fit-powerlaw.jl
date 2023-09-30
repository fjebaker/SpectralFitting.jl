using Test
using SpectralFitting

include("../dummies.jl")

# generate some fake powerlaw data
dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = "count / s keV")

# test that both julia and xspec implementations can fit simple
model = PowerLaw()
prob = FittingProblem(model, dummy_data)

# test inference
@inferred SpectralFitting._unpack_fitting_configuration(prob);

result = fit(prob, LevenbergMarquadt())

@test result.u ≈ [12.06629478087094, 3.080992500319396] atol = 1e-4
model = XS_PowerLaw()
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test result.u ≈ [12.06629478087094, 3.080992500319396] atol = 1e-4

# composite power law
dummy_data = make_dummy_dataset((E) -> (E^(-3.0) + E^(-1)); units = "count / s keV")

# julia implementation
model = PowerLaw() + PowerLaw(a = FitParam(1.0))
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
# for these the order keeps going wrong on the CI, so we'll just check the sum
@test result.u ≈ [10.374890467986033, 1.0, 12.277664914150247, 3.225587384436137] atol =
    1e-1

# xpsec implementation
model = XS_PowerLaw() + XS_PowerLaw(a = FitParam(1.0))
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test result.u ≈
      [11.152658761725311, 1.0372890868418811, 11.219122689264507, 3.310562124158435] atol =
    1e-1

# mixed implementation
model = XS_PowerLaw() + PowerLaw(a = FitParam(1.0))
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test result.u ≈
      [11.219123015018498, 3.310562096287135, 11.152658542859502, 1.0372890769565792] atol =
    1e-1
