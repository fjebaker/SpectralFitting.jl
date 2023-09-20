# Models index

Models wrapped from XSPEC implementations are prefixed with `XS_*`, whereas pure-Julia models are simply named, e.g. [`XS_PowerLaw`](@ref) in XSPEC vs [`PowerLaw`](@ref) in Julia.

The available models are
```@index
Pages = ["models.md"]
Order = [:type]
```

## Julia models

```@docs
PowerLaw
BlackBody
```

### Defining new models

To define a new Julia model, we need only implement the [`AbstractSpectralModel`](@ref) interface. As a pedagogical example, consider an implementation of [`PowerLaw`](@ref):

```julia
# define the model and parameters
@with_kw struct PowerLaw{F1,F2} <: AbstractSpectralModel
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Photon index."
    a::F2 = FitParam(0.5)
end

# tell the package what kind of model this is
modelkind(::Type{<:PowerLaw}) = Additive()

# define how the model acts: 
# all of the parameters are passed, in order, as arguments
function SpectralFitting.invoke!(flux, energy, ::Type{<:PowerLaw}, a)
    α = 1 - a
    α⁻¹ = inv(α)
    SpectralFitting.finite_diff_kernel!(flux, energy) do E
        α⁻¹ * E^α
    end
end
```

Note that each parameter has its own parametric type, here `::F1` and `::F2`. This is so that the type of fit parameters can be changed, or alternative number types inserted (see [Why & how](@ref)). The standard parameter type is [`FitParam`](@ref), which is what we have used here.

The model is defined to be additive by overloading [`modelkind`](@ref).

!!! danger
    Additive models _must have_ a normalisation parameter with the symbol `K`, which is, however, _not passed_ to [`SpectralFitting.invoke!`](@ref). Multiplicative and convolutional models have no such requirements.

We have also used the [`SpectralFitting.finite_diff_kernel!`](@ref) utility function, which is designed to help us efficiently compute finite-difference integrations of flux, just as XSPEC does.

!!! note
    To add additional XSPEC models, see [Wrapping new XSPEC models](@ref).

## XSPEC models

XSPEC models frequently have tabular data dependencies, without which the models fail to invoke (see [Model data availability](@ref)). If the data files are known but not present, the XSPEC models will throw an error with instructions for downloading the missing data. If the data files are unknown, Julia may crash catastrophically. If this is the case, often a single line will be printed with the LibXSPEC error, specifying the name of the missing source file. This can be registered as a data dependency of a model using [`SpectralFitting.register_model_data`](@ref).

The first time any XSPEC model is invoked, SpectralFitting checks to see whether requisite data is needed, and whether the data is downloaded. Subsequent calls will hit a lookup cache instead to avoid run-time costs of performing this check.

```@docs
XS_PowerLaw
XS_BlackBody
XS_BremsStrahlung
XS_Gaussian
XS_Laor
XS_DiskLine
XS_PhotoelectricAbsorption
XS_WarmAbsorption
XS_CalculateFlux
XS_KerrDisk
XS_KyrLine
```

### Wrapping new XSPEC models

SpectralFitting exports a helpful macro to facilitate wrapping additional XSPEC models.

```@docs
@xspecmodel
SpectralFitting.register_model_data
```

## Generating model fingerprints

To generate the unicode plot to add as a fingerprint, we use a simple function:

```@example
using SpectralFitting, UnicodePlots

function plotmodel(energy, model)
    flux = invokemodel(energy, model)
    lineplot(
        energy[1:end-1], 
        flux, 
        title=String(SpectralFitting.FunctionGeneration.model_base_name(typeof(model))), 
        xlabel="E (keV)", 
        canvas=DotCanvas
    )
end

# e.g. for XS_PowerLaw()
energy = collect(range(0.1, 20.0, 100))
plotmodel(energy, XS_PowerLaw())
```

