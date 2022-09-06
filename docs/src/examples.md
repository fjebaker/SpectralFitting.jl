# Spectral fitting examples

```@setup sf_examples 
using SpectralFitting
using Plots
ENV["GKSwstype"]="nul"
Plots.default(show=false)
```

Below are a number of examples illustrating what this package may be used for.

## Using the model library

The model library details a model algebra (see [`AbstractSpectralModelKind`](@ref)). An example use of this may be to construct a complex model from a series of simpler models, and invoke the models on a given energy grid:

```@example sf_examples 
using SpectralFitting
using Plots 

model = PhotoelectricAbsorption() * (PowerLaw() + BlackBody()) 

# define energy grid
energy = collect(range(0.1, 12.0, 100))

flux = invokemodel(energy, model)

plot(energy[1:end-1], flux)
```

The full model library of available models is listed in [Model index](./models.md).
