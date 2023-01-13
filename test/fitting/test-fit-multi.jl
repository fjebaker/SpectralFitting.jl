using Test
using SpectralFitting

include("../dummies.jl")

# generate some fake powerlaw data with three components
dummy_data = make_dummy_dataset((E) -> (E^(-0.1) + E^(-3.0) + E^(-1.0)))

# model with two components
model = PowerLaw() + PowerLaw()
prob = FittingProblem(MultiModel(model, model), MultiDataset(dummy_data, dummy_data))

result = fit(prob, LevenbergMarquadt())

# both models should fit more or less the same
@test sum(result.results[1].u) ≈ sum(result.results[2].u) atol = 1e-2
@test result.results[1].χ2 ≈ result.results[2].χ2 atol = 1e-2

# now change second models data to ensure the normalisations fit independelty 
dummy_data2 = deepcopy(dummy_data)
dummy_data2._data .*= 3.0

prob = FittingProblem(MultiModel(model, model), MultiDataset(dummy_data, dummy_data2))

result = fit(prob, LevenbergMarquadt())
# photon indices should still be the same
r1 = result.results[1]
r2 = result.results[2]
@test r1.u[2] + r1.u[4] ≈ r2.u[2] + r2.u[4]

# now we bind one of the normalisations
bind!(prob, :K_1)

result = fit(prob, LevenbergMarquadt())
# ensure bound parameter is identical
r1 = result.results[1]
r2 = result.results[2]
@test r1.u[1] == r2.u[1]

# now the second model will be different
model2 = PowerLaw() + PowerLaw() + PowerLaw()

prob = FittingProblem(MultiModel(model, model2), MultiDataset(dummy_data, dummy_data2))

# ensure we can fit this fine
result = fit(prob, LevenbergMarquadt())

@test result.results[1].χ2 ≈ 2.86 atol = 1e-2
@test result.results[2].χ2 ≈ 0.079 atol = 1e-3

# and we can still bind
bind!(prob, :K_1)

result = fit(prob, LevenbergMarquadt())

@test result.results[1].χ2 ≈ 97.5 atol = 1.0
@test result.results[2].χ2 ≈ 95.0 atol = 1.0
