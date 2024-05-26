# Model index

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
        title=String(Base.typename(typeof(model)).name), 
        xlabel="E (keV)", 
        canvas=DotCanvas
    )
end

# e.g. for XS_PowerLaw()
energy = collect(range(0.1, 20.0, 100))
plotmodel(energy, XS_PowerLaw())
```

