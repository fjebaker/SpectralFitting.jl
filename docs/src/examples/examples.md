# Spectral fitting examples

```@setup sf_examples 
using SpectralFitting
using Plots
ENV["GKSwstype"]="nul"
Plots.default(show=false)
```

Below are a number of examples illustrating how this package may be used.

## Using the model library

The model library details a model algebra (see [`AbstractSpectralModelKind`](@ref)) for composing models together. An example use of this may be to construct a complex model from a series of simpler models, and invoke the models on a given energy grid:

```@example sf_examples 
using SpectralFitting
using Plots 

model = PhotoelectricAbsorption() * (PowerLaw() + BlackBody()) 

# define energy grid
energy = collect(range(0.1, 12.0, 100))

flux = invokemodel(energy, model)

plot(energy[1:end-1], flux)
```

Note this energy grid may be arbitrarily spaced, but, like XSPEC, assumes the bins are contiguous, i.e. that the high energy limit of one bin is the low energy limit of the next.

The full model library of available models is listed in [Model index](./models.md).
