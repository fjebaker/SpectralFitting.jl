using Test
using SpectralFitting
import Random

lower_bounds = (0.1, 1e-3)
upper_bounds = (20.0, 3.0)

energy = collect(range(0.1, 20.0, 100))
model = XS_PhotoelectricAbsorption()

flux = similar(energy)[1:end-1]

surrogate = make_surrogate_function(
    (x, y) -> SpectralFitting.RadialBasis(x, y, lower_bounds, upper_bounds),
    model,
    lower_bounds,
    upper_bounds,
)

# number of points the surrogate has been trained on
length(surrogate.surrogate.x)

Random.seed!(1)

nh = 2.0
model.ηH.value = nh

f = invokemodel(energy, model)

f̂ = map(energy[1:end-1]) do e
    v = (e, nh)
    surrogate.surrogate(v)
end

@test sum(f .- f̂) ≈ -0.299 atol=1e-3

optimize_accuracy!(surrogate; maxiters = 50)

f̂ = map(energy[1:end-1]) do e
    v = (e, nh)
    surrogate.surrogate(v)
end

@test sum(f .- f̂) ≈ 0.196 atol=1e-3

sm = make_model(surrogate)

# smoke test these
invokemodel(collect(range(0.01, 10.0, 100)), sm)
new_model = sm * PowerLaw()
invokemodel(collect(range(0.01, 10.0, 1000)), new_model)
