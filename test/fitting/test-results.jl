using Test, SpectralFitting

include("../dummies.jl")

dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
model = PowerLaw()
prob = FittingProblem(model, dummy_data)

# test inference
@inferred SpectralFitting._unpack_fitting_configuration(prob);

result = fit(prob, LevenbergMarquadt())

@test result[1].u ≈ [12.066, 3.0810] atol = 1e-3

@test measure(ChiSquared(), result) ≈ result.χ2 rtol = 1e-3
@test measure(ChiSquared(), result, [11.0, 2.0]) ≈ 281423.5 rtol = 1e-3
# can measure different parameters
@test measure(ChiSquared(), result[1], [11.0, 2.0]) ≈ 281423.5 rtol = 1e-3

# check these don't error
@test invoke_result(result) ≈ invoke_result(result[1])
invoke_result(result, [11.0, 2.0])

# mismatch
@test_throws "" invoke_result(result, [11.0])
