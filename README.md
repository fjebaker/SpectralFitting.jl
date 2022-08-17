<p align="center">
  <img width="50%" alt="BRImage" src="docs/src/assets/logo.svg">
</p>

# SpectralFitting.jl

<a href="https://fjebaker.github.io/SpectralFitting.jl/dev/">
<img alt="Docs" src="https://img.shields.io/badge/docs-dev-blue.svg"/>
</a>

<p align="center"> <i> This package is in development and subject to rapid breaking changes.</i> </p>

## Setup

Install LibXSPEC_jll:
```julia
julia>] add "https://github.com/astro-group-bristol/LibXSPEC_jll.jl#master"
```

This package contains pre-built binaries of the XSPEC models, along with the data-files, and is consequently quite big. A future build will make it available for more Julia versions, more OSs, and remove the data-files, such that they are downloaded _as-needed_.

Next, install this package:
```julia
julia>] add "https://github.com/fjebaker/SpectralFitting.jl"
```

Check everything went okay by importing the package in the REPL:
```julia
julia> using SpectralFitting
# Solar Abundance Vector set to angr:  Anders E. & Grevesse N. Geochimica et Cosmochimica Acta 53, 197 (1989)
# Cross Section Table set to vern:  Verner, Ferland, Korista, and Yakovlev 1996
```

## Using the models

Julia types which wrap XSPEC implementations of models are prefixed with `XS_*`, whereas the pure-Julia alternatives will just be the model name, e.g. `XS_PowerLaw` vs `PowerLaw`.

Here is an example of a composite model, equivalent to `phabs(powerlaw)` in XSPEC:
```julia
using SpectralFitting
using Plots

model = XS_PhotoelectricAbsorption() * XS_PowerLaw()

energy = collect(range(0.1, 20.0, 200))
fluxes = make_fluxes(energy, flux_count(model))

# compiles specialized function and calls XSPEC models
invokemodel!(fluxes, energy, model)

outflux = first(fluxes)

# energy is an array of bin edges, which contains one more element than flux
plot(energy[1:end-1], outflux)
```

## Defining new models

To wrap a different XSPEC model, SpectralFitting includes a handy macro, with the following syntax:
```julia
@xspecmodel MODEL_KIND C_FUNCTION struct MODEL_NAME{F1,F2,...}
    param1::F1 = FitParam(1.0)
    param2::F2 = FrozenFitParam(1.0)
    # ...
end
```

- `MODEL_KIND` may be `Additive`, `Multiplicative`, or `Convolutional`. In the case of `Additive`, models must specify a normalisation constant `K` as their **first parameter**.
- `C_FUNCTION` is a symbol which describes the corresponding XSPEC function to call.

An example:

```julia
SpectralFitting.@xspecmodel Additive :C_bbody struct XS_BlackBody{F1,F2}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Temperature (keV)."
    T::F2 = FitParam(3.0)
end
```

The macro will generate all the necessary `ccall` wrapping and additional model definitions for you.

**Note**: not every model will work out-of-the-box, due to an issue with Fortran versions. LibXSPEC_jll has only been compiled against Fortran 5.x, whereas the Bristol servers run either Fortran 4.x or Fortran >7.x. This will be patched soon, to ship the models for a wider suite of compiler versions.

To create a Julia model, we define a similar structure, but specify the model kind and action separately, e.g.:

```julia
@with_kw struct BlackBody{F1,F2} <: AbstractSpectralModel
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Temperature (keV)."
    kT::F2 = FitParam(2.0)
end

# specify model kind as a trait
modelkind(::Type{<:BlackBody}) = Additive()

# define action in-place
@fastmath function SpectralFitting.invoke!(flux, energy, ::Type{<:BlackBody}, kT)
    finite_diff_kernel!(flux, energy) do E
        8.0525 * E^2 / (kT^4 * (exp(E / kT) - 1))
    end
end
```

## Other currently undocumented features

- MCMC fitting
- Speedy surrogate models
- Fits using automatic differentiation