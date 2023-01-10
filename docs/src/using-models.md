# Using spectral models

```@setup using_models
using SpectralFitting
using Plots
ENV["GKSwstype"]="nul"
Plots.default(show=false)
```

Models exist as three different kinds, defined by an [`AbstractSpectralModelKind`](@ref) trait.
```@docs
AbstractSpectralModelKind
Additive
Multiplicative
Convolutional
```

## Model abstraction

All spectral models are a sub-type of [`AbstractSpectralModel`](@ref), including [`CompositeModel`](@ref) and [`SurrogateSpectralModel`](@ref).

```@docs
AbstractSpectralModel
SpectralFitting.invoke!
modelkind
```

## Model data availability

Many of the XSPEC implemented models use tabular data, such as FITS, and return results interpolated from these pre-calculated tables. In some cases, these table models have data files that are multiple gigabytes in size, and would be very unwieldy to ship indiscriminantly. SpectralFitting attempts to circumnavigate this bloat by downloading the model data on an _ut opus_ basis.

```@docs
SpectralFitting.download_model_data
SpectralFitting.download_all_model_data
```

Special care must be taken if new XSPEC models are wrapped to ensure the data is available. For more on this, see [Wrapping new XSPEC models](@ref).

Model data may also alternatively be copied in _by-hand_ from a HEASoft XSPEC source directory. In this case, the location to copy the data to may be determined via `joinpath(SpectralFitting.LibXSPEC_jll.artifact_dir, "spectral", "modelData")`.

## Instantiating and invoking models

Models may be composed together to create more complex spectra, with the algebra defined by the [`AbstractSpectralModelKind`](@ref). 

For example, a [`XS_PowerLaw`](@ref) model may multiply element-wise by a [`XS_PhotoelectricAbsorption`](@ref) model:

```@example using_models
model = XS_PhotoelectricAbsorption() * XS_PowerLaw()
model
```

A model may be evaluated by using [`invokemodel`](@ref):

```@example using_models
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, model)
```

!!! note
    XSPEC models are implemented to accept a single energy vector with all of the low and high energy bin edges, and return a flux array with the flux in each energy bin. For consistency, SpectralFitting implement models in the same way, and therefore it is always the case that:
    ```julia
    length(flux) == length(energy) - 1
    ```

```@docs
invokemodel
invokemodel!
```

## Utilities

SpectralFitting defines a number of utility functions for model introspection and manipulation.

```@docs
flux_count
get_param
get_param_count
get_param_types
get_param_symbols
get_param_symbol_pairs
get_params
get_params_value
```