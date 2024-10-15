using Test
using SpectralFitting
using OptimizationOptimJL

include("../dummies.jl")

# generate some fake powerlaw data with three components
dummy_data = make_dummy_dataset(
    (E) -> (E^(-0.1) + E^(-3.0) + E^(-1.0));
    units = u"counts / (s * keV)",
)

# model with two components
model = PowerLaw(K = FitParam(10.0)) + PowerLaw(K = FitParam(10.0))
prob = FittingProblem(
    FittableMultiModel(model, model),
    FittableMultiDataset(dummy_data, dummy_data),
)

# check inferred types for multi models
@inferred SpectralFitting._build_parameter_mapping(prob.model, prob.bindings)
@inferred SpectralFitting._build_parameter_mapping(prob.model, prob.bindings)
@inferred SpectralFitting._unpack_config(prob)

result = fit(prob, NelderMead(), autodiff = SpectralFitting.Optimization.SciMLBase.NoAD())

# both models should fit more or less the same
@test sum(result[1].u) ≈ sum(result[2].u) atol = 1e-2
@test result[1].χ2 ≈ result[2].χ2 atol = 1e-2
@test sum(result.χ2s) ≈ 5.7238519700611725 atol = 1e-2

# now with the AD backend
result = fit(prob, BFGS())

# both models should fit more or less the same
@test sum(result[1].u) ≈ sum(result[2].u) atol = 1e-2
@test result[1].χ2 ≈ result[2].χ2 atol = 1e-2
@test sum(result.χ2s) ≈ 5.7238519700611725 atol = 1e-2

# do a single model

prob = FittingProblem(model => dummy_data)

result = fit(prob, BFGS())

# both models should fit more or less the same
@test result.u ≈ [14.970, 3.0974, 17.254, 0.285] atol = 1e-2
@test result.χ2 ≈ 2.8619 atol = 1e-2

# now with different statistic

prob = FittingProblem(model => dummy_data)
result = fit(prob, BFGS(); stat = Cash())

# both models should fit more or less the same
@test result.u ≈ [12.217, 3.288, 19.053, 0.339] atol = 1e-2
@test result.χ2 ≈ 0.305 atol = 1e-2
