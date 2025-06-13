using Test
using SpectralFitting

include("../dummies.jl")

dataset = make_dummy_dataset(i -> i^-2)

struct TestModelWrapper{M,T,K} <: SpectralFitting.AbstractModelWrapper{M,T,K}
    model::M
end

TestModelWrapper(model::AbstractSpectralModel{T,K}) where {T,K} =
    TestModelWrapper{typeof(model),T,K}(model)

# single wrapper
model = TestModelWrapper(PowerLaw())
domain = collect(range(0.1, 10.0, 100))
@inferred invokemodel(domain, model)

# composite wrapper
model = TestModelWrapper(PowerLaw() + PowerLaw())
domain = collect(range(0.1, 10.0, 100))

output = @inferred invokemodel(domain, model)
@test size(allocate_model_output(model, domain)) == (99, 2)
@test size(SpectralFitting.parameter_vector(model)) == (4,)

backing_output = invokemodel(domain, SpectralFitting.backing_model(model))

@test output ≈ backing_output rtol = 0.01

@inferred SpectralFitting.parameter_names(typeof(model))

prob = FittingProblem(model => dataset)

conf = @inferred FittingConfig(prob)

# some smoke tests
@inferred SpectralFitting.calculate_objective!(conf, [1.0, 1.0, 1.0, 1.0])

fit(prob, LevenbergMarquadt())

# composite wrapper but as part of another composite model
model = TestModelWrapper(PowerLaw() + PowerLaw()) + PowerLaw()
domain = collect(range(0.1, 10.0, 100))
@inferred invokemodel(domain, model)

prob = FittingProblem(model => dataset)

fit(prob, LevenbergMarquadt())
