# Surrogate models

```@setup surrogate_example
using SpectralFitting
using Plots
ENV["GKSwstype"]="nul"
Plots.default(show=false)
```

Surrogate models allow you to create fast or memory efficient approximations of model components, or assist in optimizing some objective function directly. SpectralFitting uses the [Surrogates.jl](https://github.com/SciML/Surrogates.jl) library of models, that can yields pure-Julia surrogate models. Consequently surrogate models also permit use of automatic differentiation in fitting, and are therefore powerful tools for improving fitting performance.

## Surrogates overview

Any function may be wrapped as a surrogate model using the [`SurrogateSpectralModel`](@ref) type.

```@docs
SurrogateSpectralModel
```

To facilitate easy surrogate builds, SpectralFitting exports a number of utility functions.

```@docs
make_surrogate_function
optimize_accuracy!
```

## Creating a surrogate for `XS_PhotoelectricAbsorption`

Before we start, let us discuss a number of benefits the use of surrogate models may bring us:

- [`SurrogateSpectralModel`](@ref) permit use of automatic differentiation.
- Surrogate models may be allocation-free depending on setup, whereas XSPEC wrappers will always have to allocate for type-conversions.
- Surrogate models may be considerably faster, especially for table models.
- Surrogate models are shareable (see [Sharing surrogate models](@ref)), and are tunable in size.

[`XS_PhotoelectricAbsorption`](@ref) is an XSPEC model that is wrapped by a thin C-wrapper into Julia. The implementation of this model is a number of Fortran routines from the late 90s, including a tabulation of ~3000 lines of data that has been copied directly into the Fortran source code.

The performance of this model represents its complexity:

```@example surrogate_example
using BenchmarkTools

model = XS_PhotoelectricAbsorption()
energy = collect(range(0.1, 20.0, 100))
flux = make_flux(energy)

@benchmark invokemodel!($flux, $energy, $model)
```

The surrogate we'll construct will have to be tailored a little to the data we wish to fit, as we need to specify the parameter ranges our surrogate should learn. For example, we might be interested in energies between ``0.1`` and ``20`` keV, with equivalent hydrogen column ``\eta H`` anywhere between ``10^{-3}`` and ``30``. We specify these bounds using tuples

```@example surrogate_example
lower_bounds = (0.1, 1e-3)
upper_bounds = (20.0, 30.0)
nothing # hide
```

!!! note
    The first index is always the energy bounds, and the subsequent indices are the parameters in the same order they are defined in the model structure.

Next, we use [`make_surrogate_function`](@ref) to build and optimize a surrogate function for our model. By default, the surrogate uses linear radial basis functions, and seeds the coefficients with a number of seed points. This function then improves the accuracy of the model using [`optimize_accuracy!`](@ref), until a maximal number of iterations has been reached.

For illustration purposes, we'll omit the accuracy improving step, and perform this ourselves. We can do this by setting `optimization_samples = 0` in the keyword arguments:

```@example surrogate_example
surrogate = make_surrogate_function(
    model, 
    lower_bounds, 
    upper_bounds
    ;
    optimization_samples = 0
)

# number of points the surrogate has been trained on
length(surrogate.x)
```

We can examine how well our surrogate reconstructs the model for a given test parameter:

```@example surrogate_example
import Random # hide
Random.seed!(1) # hide
# random test value
ηh_test = rand(1.0:30.0)

f = invokemodel(energy, model, ηh_test)

f̂ = map(energy[1:end-1]) do e
    v = (e, ηh_test)
    surrogate(v)
end
p = plot(energy[1:end-1], f, label="model", legend=:bottomright, xlabel="E (keV)") # hide
plot!(energy[1:end-1], f̂, label="surr") # hide
p # hide
```

Now we'll use [`optimize_accuracy!`](@ref) to improve the faithfulness of our surrogate. This requires making use of [`wrap_model_as_objective`](@ref) as a little wrapper around our model:

```@example surrogate_example
obj = wrap_model_as_objective(model)

optimize_accuracy!(surrogate, obj, lower_bounds, upper_bounds; maxiters=250)

length(surrogate.x)
```

We can plot the surrogate model again and see the improvement.
```@example surrogate_example
f̂2 = map(energy[2:end]) do e
    v = (e, ηh_test)
    surrogate(v)
end
plot!(energy[1:end-1], f̂2, label="surr+") # hide
p # hide
```

Tight. We can also inspect the memory footprint of our model:

```@example surrogate_example
# in bytes
Base.summarysize(surrogate)
```
This may be reduced by lowering `maxiters` in [`optimize_accuracy!`](@ref) at the cost of decreasing faithfulness. However, compare this to the Fortran tabulated source file in the XSPEC source code, which is approximately 224 Kb -- about 22x larger. The surrogate models are considerably more portable at this level.

Inspecting the plot shows that there is a small domain near zero energy where the surrogate is calculating a negative value. We can either try to improve the accuracy to fix this, or specify a clamping function: for [`XS_PhotoelectricAbsorption`](@ref), we know the multiplicative factor is constrained between 0 and 1. Therefore, we may adjust our surrogate to reflect this:

```@example surrogate_example
clamped_surrogate = (v) -> clamp(surrogate(v), 0.0, 1.0)
nothing # hide
```

## Using a surrogate spectral model

Now that we have the surrogate model, we use [`SurrogateSpectralModel`](@ref) to wrap it into an [`AbstractSpectralModel`](@ref). The constructor also needs to know the model kind, have a copy of the model parameters, and know which symbols to represent the parameters with.

```@example surrogate_example
sm = SurrogateSpectralModel(
    Multiplicative(),
    clamped_surrogate,
    (FitParam(1.0),), # must be a tuple
    (:ηh,) # must also be a tuple
)
```

We can now use the familiar API and attempt to benchmark the performance:

```@example surrogate_example
@benchmark invokemodel!($flux, $energy, $sm)
```

These allocations are coming from the global closure we have in the  `clamped_surrogate` lambda. We can actually elide these using a proper closure:

```@example surrogate_example
make_surr(surrogate) = SpectralFitting.SurrogateSpectralModel(
    Multiplicative(),
    (v) -> clamp(surrogate(v), 0.0, 1.0), 
    (FitParam(1.0),), 
    (:ηH,)
)

sm2 = make_surr(surrogate)
@benchmark invokemodel!($flux, $energy, $sm2)
```

Comparing this to the initial benchmark of [`XS_PhotoelectricAbsorption`](@ref), we see about a 60x speedup, with no allocations. Furthermore, this surrogate model is now automatic differentiation ready!

## Sharing surrogate models

To export and import surrogate models, [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) is recommended.