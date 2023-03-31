using Test
using SpectralFitting

struct DummySimplePowerLaw{T,F} <: AbstractSpectralModel{T,Additive}
    K::T
    m::T
end
function DummySimplePowerLaw(; K = FitParam(1.0), m = FitParam(1.0))
    DummySimplePowerLaw{typeof(K),SpectralFitting.FreeParameters{(:K, :m)}}(K, m)
end
function SpectralFitting.invoke!(out, x, model::DummySimplePowerLaw)
    @. out = x + model.m
end
SpectralFitting.Δoutput_length(::Type{<:DummySimplePowerLaw}) = 0

# generate a powerlaw simple data
x = collect(range(1.0, 10.0, 18))
y = @. 13.12 * x + 42
data = SimpleDataset("example", x, y)

model = DummySimplePowerLaw()

prob = FittingProblem(model, data)

res = fit(prob, LevenbergMarquadt())
@test res.u[1] ≈ 13.12 atol = 0.1
@test res.u[2] ≈ 3.201 atol = 0.1

# test error constructors
data = SimpleDataset("example", x, y, x_err = 0.1 .* x, y_err = 0.1 * y)
