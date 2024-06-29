using SpectralFitting
using Test

include("../dummies.jl")

# put a couple of delta emission lines together
lines = DeltaLine(; E = FitParam(3.0), K = FitParam(2.0)) + DeltaLine(; E = FitParam(7.0))

# construct the convolutional wrapper
base_model = GaussianLine(; μ = FitParam(1.0), σ = FitParam(0.3))
conv = AsConvolution(base_model)

model = conv(lines)

domain = collect(range(0.0, 10.0, 150))

plot(domain[1:end-1], invokemodel(domain, lines))
plot(domain[1:end-1], invokemodel(domain, model))

output = invokemodel(domain, model)

@test sum(output) ≈ 3.2570820013702395 atol = 1e-4
@test output[10] ≈ 0.0036345342427057687 atol = 1e-4
@test output[40] ≈ 0.055218163108951814 atol = 1e-4

# simulate a model spectrum
dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
sim = simulate(model, dummy_data; seed = 42)

model.μ_1.frozen = true
model.K_1.frozen = true
model.K_2.frozen = true
model.E_1.frozen = true
model.E_2.frozen = true

# change the width
model.σ_1.value = 0.1
model

begin
    prob = FittingProblem(model => sim)
    result = fit(prob, LevenbergMarquadt())
end
@test result.χ2 ≈ 76.15221077389369 atol = 1e-3

# put a couple of delta emission lines together
lines =
    DeltaLine(; E = FitParam(3.0), K = FitParam(2.0), width = 0.1) +
    DeltaLine(; E = FitParam(7.0))
model = conv(lines)

sim = simulate(model, dummy_data; seed = 42)

# now see if we can fit the delta line
model.μ_1.frozen = true
model.K_1.frozen = true
model.K_2.frozen = true
model.E_1.frozen = true
model.E_2.frozen = true
model.σ_1.frozen = true

model.E_2.frozen = false
model.E_2.value = 2.0
model.K_2.frozen = true
# model.K_2.value = 2.0

model
begin
    prob = FittingProblem(model => sim)
    result = fit(prob, LevenbergMarquadt(); verbose = true)
end
@test result.χ2 ≈ 75.736 atol = 1e-3
