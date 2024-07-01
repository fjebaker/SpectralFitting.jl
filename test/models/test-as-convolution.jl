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

output = invokemodel(domain, model)

@test sum(output) ≈ 4.050608829695485 atol = 1e-4
@test output[10] ≈ 0.0022170544439135222 atol = 1e-4
@test output[40] ≈ 0.058630601782812125 atol = 1e-4

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
@test result.χ2 ≈ 76.71272868245076 atol = 1e-3
@test result.u[1] ≈ 0.3 atol = 1e-2

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
@test result.χ2 ≈ 76.66970981760741 atol = 1e-3
@test result.u[1] ≈ 3.0 atol = 1e-2
