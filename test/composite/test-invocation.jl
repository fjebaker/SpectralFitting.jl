using Test
using SpectralFitting

include("../dummies.jl")

model = DummyAdditive()

# just check that it has no errors
checker_free_frozen(model) = SpectralFitting.FunctionGeneration.generated_model_call!(
    Vector{Float64},
    Vector{Float64},
    typeof(model),
    Vector{Float64},
    Vector{Float64},
)
checker_all(model) = SpectralFitting.FunctionGeneration.generated_model_call!(
    Vector{Float64},
    Vector{Float64},
    typeof(model),
    Vector{Float64},
)
_ = checker_free_frozen(model)
_ = checker_free_frozen(model + model)
_ = checker_all(model)
_ = checker_all(model + model)

# check that the generated functions work too

energy = collect(range(0.1, 100.0, 100))
# single model
out_flux = invokemodel(energy, model)
@test all(out_flux .== 6)

# added models
out_flux = invokemodel(energy, model + model)
@test all(out_flux .== 12)

# big added models
out_flux = invokemodel(
    energy,
    model +
    model +
    model +
    model +
    model +
    model +
    model +
    model +
    model +
    model +
    model +
    model +
    model +
    model +
    model +
    model,
)
@test all(out_flux .== 96)

# test multiplicative with additive
cm = DummyMultiplicative() * DummyAdditive()
out_flux = invokemodel(energy, cm)
@test all(out_flux .== 30)

# test inplace version
invokemodel!(out_flux, energy, DummyAdditive())
@test all(out_flux .== 6)