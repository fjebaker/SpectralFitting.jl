using Test
using SpectralFitting
using OptimizationOptimJL

include("../dummies.jl")

# generate some fake powerlaw data with three components
dummy_data = make_dummy_dataset((E) -> (E^(-0.1) + E^(-3.0) + E^(-1.0)))

# model with two components
model = PowerLaw() + PowerLaw()
prob = FittingProblem(MultiModel(model, model), MultiDataset(dummy_data, dummy_data))

result = fit(prob, ChiSquared(), NelderMead())
#
# both models should fit more or less the same
@test sum(result.results[1].u) ≈ sum(result.results[2].u) atol = 1e-2
@test result.results[1].χ2 ≈ result.results[2].χ2 atol = 1e-2

