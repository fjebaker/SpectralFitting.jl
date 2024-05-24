# Surrogate models

Surrogate models allow you to create fast or memory efficient approximations of model components, or assist in optimizing some objective function directly. SpectralFitting uses the [Surrogates.jl](https://github.com/SciML/Surrogates.jl) library of models, that yields pure-Julia surrogate models. Consequently, surrogate models also permit use of automatic differentiation in fitting, and are therefore powerful tools for improving fitting performance.

## Surrogates overview

!!! warning

    The surrogate model optimization does not work well for most XSPEC models currently. This is being actively developed.

Any function may be wrapped as a surrogate model using the [`SurrogateSpectralModel`](@ref) type.

```@docs
SurrogateSpectralModel
```

To facilitate easy surrogate builds, SpectralFitting exports a number of utility functions.

```@docs
make_surrogate_harness
optimize_accuracy!
```

## Creating a surrogate for `XS_PhotoelectricAbsorption`

Before we start, let us discuss a number of benefits the use of surrogate models may bring us:

- [`SurrogateSpectralModel`](@ref) permit use of automatic differentiation.
- Surrogate models may be allocation-free depending on setup, whereas XSPEC wrappers will always have to allocate for type-conversions.
- Surrogate models may be considerably faster, especially for table models.
- Surrogate models are shareable (see [Sharing surrogate models](@ref)), and are tunable in size.

[`XS_PhotoelectricAbsorption`](@ref) is an XSPEC model that is wrapped by a thin C-wrapper into Julia. The implementation of this model is a number of Fortran routines from the late 90s, including a tabulation of ~3000 lines of data that has been copied directly into the Fortran source code.

The performance of this model represents its complexity.

```@example surrogate_example
using SpectralFitting 

energy = collect(range(0.1, 20.0, 200))
model = XS_PhotoelectricAbsorption()

flux = similar(energy)[1:end-1]
```

Benchmarking with [BenchmarkTools.jl](https://juliaci.github.io/BenchmarkTools.jl/stable/):

```@example surrogate_example
using BenchmarkTools
@btime invokemodel!($flux, $energy, $model);
```

The surrogate we'll construct will have to be tailored a little to the data we wish to fit, as we need to specify the parameter ranges our surrogate should learn. For example, we might be interested in energies between ``0.1`` and ``20`` keV (expressed in our domain), with equivalent hydrogen column ``\eta``H anywhere between ``10^{-3}`` and ``30``. We specify the parameter bounds using tuples:

```@example surrogate_example
lower_bounds = (1e-3,)
upper_bounds = (30.0,)
nothing # hide
```

!!! note
    The first index is always the energy bounds, and the subsequent indices are the parameters in the same order they are defined in the model structure.

Next, we use [`make_surrogate_harness`](@ref) to build and optimize a surrogate function for our model. By default, the surrogate uses linear radial basis functions, and seeds the coefficients with a number of seed points. This function then improves the accuracy of the model using [`optimize_accuracy!`](@ref), until a maximal number of iterations has been reached.

For illustration purposes, we'll omit the accuracy improving step, and perform this ourselves. We can do this by setting `optimization_samples = 0` in the keyword arguments:

```@example surrogate_example
using Surrogates

harness = make_surrogate_harness(
    (x, y) -> RadialBasis(x, y, lower_bounds, upper_bounds),
    energy,
    model,
    lower_bounds,
    upper_bounds,
)

# number of points the surrogate has been trained on
length(harness.surrogate.x)
```

We can examine how well our surrogate reconstructs the model for a given test parameter:

```@example surrogate_example
using Plots
import Random # hide
Random.seed!(1) # hide
# random test value
ηh_test = rand() * 3 + lower_bounds[1]

model.ηH.value = ηh_test
f = invokemodel(energy, model)

f̂ = harness.surrogate([ηh_test])

p = plot(energy[1:end-1], f, label="model", legend=:bottomright, xlabel="E (keV)") # hide
plot!(energy[1:end-1], f̂, label="surr") # hide
p # hide
```

Now we'll use [`optimize_accuracy!`](@ref) to improve the faithfulness of our surrogate. This requires making use of [`wrap_model_as_objective`](@ref) as a little wrapper around our model:

```@example surrogate_example
optimize_accuracy!(harness; maxiters=50)

length(harness.surrogate.x)
```

We can plot the surrogate model again and see the improvement.
```@example surrogate_example
new_f̂ = harness.surrogate([ηh_test])
plot!(energy[1:end-1], new_f̂, label="surr+") # hide
p # hide
```

Tight. We can also inspect the memory footprint of our model:

```@example surrogate_example
# in bytes
Base.summarysize(harness) 
```
This may be reduced by lowering `maxiters` in [`optimize_accuracy!`](@ref) at the cost of decreasing faithfulness. However, compare this to the Fortran tabulated source file in the XSPEC source code, which is approximately 224 Kb. The surrogate model with all it's training data is of the same order.

## Using a surrogate spectral model

Now that we have the surrogate model, we use [`SurrogateSpectralModel`](@ref) to wrap it into an [`AbstractSpectralModel`](@ref). The constructor also needs to know the model kind, have a copy of the model parameters, and know which symbols to represent the parameters with.

```@example surrogate_example
sm = make_model(harness)
```

We can now use the familiar API and attempt to benchmark the performance:

```@example surrogate_example
@btime invokemodel!($flux, $energy, $sm);
```

Comparing this to the initial benchmark of [`XS_PhotoelectricAbsorption`](@ref), we see about a significant speedup, with no allocations, _and_ this surrogate model is now automatic differentiation ready.

## Evaluating the model

```@example surrogate_example
p_range = collect(range(1.0, 30.0))

fluxes_vecs = map(p_range) do p
    model.ηH.value = p
    f = invokemodel(energy, model)
end
fluxes_mat = reduce(hcat, fluxes_vecs)

surface(p_range, energy[1:end-1], fluxes_mat, xlabel = "ηH", ylabel = "E", zlabel = "f", title = "Model")
```

```@example surrogate_example
s_fluxes_vecs = map(p_range) do p
    sm.params[1].value = p
    display(sm)
    f = invokemodel(energy, sm)
end
s_fluxes_mat = reduce(hcat, s_fluxes_vecs)

surface(p_range, energy[1:end-1], s_fluxes_mat, xlabel = "ηH", ylabel = "E", zlabel = "f", title = "Surrogate")
```

## Sharing surrogate models

To export and import surrogate models, [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) is recommended.