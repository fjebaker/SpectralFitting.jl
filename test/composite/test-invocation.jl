using Test
using SpectralFitting

include("../dummies.jl")
include("../fuzz.jl")

model = DummyAdditive()

# test single invocation
energy = collect(range(0.1, 100.0, 100))
flux = zeros(Float64, length(energy) - 1)

out_flux = invokemodel!(flux, energy, model)
@test all(out_flux .== 6)

# make sure normalisation is handled correctly
model2 = DummyAdditive(K = FitParam(2.0))

out_flux = invokemodel!(flux, energy, model2)
@test all(out_flux .== 12)

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

# some other special cases
_ = checker_all(DummyMultiplicative() * DummyAdditive())
_ = checker_free_frozen(DummyMultiplicative() * DummyAdditive())

for model in FUZZ_ALL_MODELS
    _ = checker_all(model)
    _ = checker_free_frozen(model)
end

# check that the generated functions work too

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
model = DummyAdditive()
invokemodel!(out_flux, energy, model)
@test all(out_flux .== 6)

# should be unchanged if called multiple times
invokemodel!(out_flux, energy, model)
invokemodel!(out_flux, energy, model)
@test all(out_flux .== 6)

# invocation with different parameters for single model
model = DummyAdditive()
free_params = [2.0, 2.0]
out_flux = invokemodel(energy, model, free_params)
@test all(out_flux .== 14.0)

flux = zeros(Float64, length(energy) - 1)
invokemodel!(flux, energy, model, free_params)
@test all(flux .== 14.0)

# multiple calls should give the same result
invokemodel!(flux, energy, model, free_params)
invokemodel!(flux, energy, model, free_params)
@test all(flux .== 14.0)

# invocation with different parameters for composite model
cm = DummyMultiplicative() * DummyAdditive()
free_params = [2.0, 2.0, 2.0]
out_flux = invokemodel(energy, cm, free_params)
@test all(out_flux .== 140.0)

flux = zeros(Float64, length(energy) - 1)
fluxes = (flux, deepcopy(flux))
invokemodel!(fluxes, energy, cm, free_params)
@test all(flux .== 140.0)

invokemodel!(fluxes, energy, cm, free_params)
invokemodel!(fluxes, energy, cm, free_params)
@test all(flux .== 140.0)

# ensure we can pass different input parameter types and have the system update accordingly
model = DummyAdditive()
free_params = Float32[2.0, 1.0]
flux = invokemodel(energy, model, free_params)
@test all(flux .== 12.0)
@test eltype(flux) == Float32

# composite
cm = DummyMultiplicative() * DummyAdditive()
free_params = Float32[2.0, 2.0, 2.0]
flux = invokemodel(energy, cm, free_params)
@test all(flux .== 140.0)
@test eltype(flux) == Float32

# composite with many frozen
cm = DummyMultiplicative() * (DummyAdditive() + DummyAdditiveWithManyFrozen())
flux = invokemodel(energy, cm)
@test all(flux .== 160.0)

# inplace variant
flux = zeros(Float64, length(energy) - 1)
fluxes = (flux, deepcopy(flux))
invokemodel!(fluxes, energy, cm)
@test all(flux .== 160.0)

# inplace with new free parameters
flux .= 0

free_params = [1.0, 2.0, 1.0, 1.0, 1.0]
invokemodel!(fluxes, energy, cm, free_params)
@test all(flux .== 160.0)
