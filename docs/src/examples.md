# Spectral fitting examples

Below are a number of examples illustrating what this package may be used for.

## Using the model library

The model library details an model algebra (see [`AbstractSpectralModelKind`](@ref)). An example use of this may be to construct a complex model from a series of simpler models, and invoke the models on a given energy grid:

```@example ex_models
using SpectralFitting
using Plots 

model = PhotoelectricAbsorption() * (PowerLaw() + BlackBody()) 

# define energy grid
energy = collect(range(0.1, 12.0, 100))

flux = invokemodel(energy, model)

plot(energy[1:end-1], flux)
```

The full model library of available models is listed in [Model index](./models.md).
