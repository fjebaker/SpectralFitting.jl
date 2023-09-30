using Test
using SpectralFitting
using OptimizationOptimJL

include("../dummies.jl")

# generate some fake powerlaw data with three components
dummy_data =
    make_dummy_dataset((E) -> (E^(-0.1) + E^(-3.0) + E^(-1.0)); units = "count / s")

# model with two components
model = PowerLaw(K = FitParam(10.0)) + PowerLaw(K = FitParam(10.0))
prob = FittingProblem(
    FittableMultiModel(model, model),
    FittableMultiDataset(dummy_data, dummy_data),
)

# todo: something is currently broken with the AD and i don't understand what
result = fit(
    prob,
    ChiSquared(),
    NelderMead(),
    autodiff = SpectralFitting.Optimization.SciMLBase.NoAD(),
)

# both models should fit more or less the same
@test sum(result[1].u) ≈ sum(result[2].u) atol = 1e-2
@test result[1].χ2 ≈ result[2].χ2 atol = 1e-2
@test sum(result.χ2s) ≈ 5.7238519700611725 atol = 1e-2
