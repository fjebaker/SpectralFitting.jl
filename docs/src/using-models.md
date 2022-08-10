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

All spectral models are a sub-type of [`AbstractSpectralModel`](@ref), including [`CompositeSpectralModel`](@ref) and [`SurrogateSpectralModel`](@ref).

```@docs
AbstractSpectralModel
SpectralFitting.invoke!
modelkind
```

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
get_all_model_params
get_all_model_params_by_value
```