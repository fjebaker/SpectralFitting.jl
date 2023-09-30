using SpectralFitting, Test

struct DummySimpleLinear{T,F} <: AbstractSpectralModel{T,Additive}
    K::T
    c::T
end
function DummySimpleLinear(; K = FitParam(1.0), c = FitParam(1.0))
    DummySimpleLinear{typeof(K),SpectralFitting.FreeParameters{(:K, :c)}}(K, c)
end
function SpectralFitting.invoke!(out, x, model::DummySimpleLinear)
    @. out = x + (model.c / model.K)
end
SpectralFitting.supports_contiguosly_binned(::Type{<:DummySimpleLinear}) = false
SpectralFitting.supports_one_to_one(::Type{<:DummySimpleLinear}) = true

x = collect(range(1.0, 10.0, 101))
y = @. 13.12 * x + 102

# first we try with one to one support
data = InjectiveData(x, y)
model = DummySimpleLinear()

@test SpectralFitting.preferred_support(model) isa SpectralFitting.OneToOne
@test size(SpectralFitting.construct_objective_cache(model, x)) == (101, 1)

@test SpectralFitting.common_support(model, data) isa SpectralFitting.OneToOne

# no erors
invokemodel(x, model)

prob = FittingProblem(model => data)

result = fit(prob, LevenbergMarquadt())
@test result.u ≈ [13.12, 102.0]

# now we try something with contiguous bins
y = @. x^-6
data = InjectiveData(x, y)
model = PowerLaw()

@test SpectralFitting.common_support(model, data) isa SpectralFitting.ContiguouslyBinned

prob = FittingProblem(model => data)

result = fit(prob, LevenbergMarquadt())
@test result.u[2] ≈ 6.1 atol = 0.1
