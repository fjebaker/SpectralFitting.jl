# Parameters

## Parameter binding

When performing a fit, it is desireable to **bind** certain parameters together. This ensures that they will have the same value; for example, if you were fitting two simultaneous datasets with two [`PowerLaw`](@ref) models, you may want to have different normalisations of the model components, but enforce the power law index to be the same. To achieve this, SpectralFitting has the [`bind!`](@ref) function that applies to your [`FittingProblem`](@ref).

```@docs
bind!
```

!!! note
    Bindings are treated not as specific to the model but specific to the [`FittingProblem`](@ref). This is because you may want to use the same model for multiple different datasets, and have slightly different binding requirements for each one (e.g. depending on the instruments you are using). If you do need the same binding applied to two different problems, you can do that with
    ```julia
    append!(prob1.bindings, prob2.bindings)
    ```
    Caution however, this will only make sense if you are using precisely the same model in both problems.


Let's try it out. We'll generate some arbitrary powerlaw spectra with different normalisations and fit them simultaneously.
```julia
using SpectralFitting, Plots

energy = collect(range(0.1, 10.0, 100))

# two different models with different normalisations
model1 = PowerLaw(K = FitParam(100.0), a = FitParam(1.2))
model2 = PowerLaw(K = FitParam(300.0), a = FitParam(1.22))

data1 = simulate(energy, model1, var = 1e-3)
data2 = simulate(energy, model2, var = 1e-3)

plot(data1, xscale = :log10, yscale = :log10)
plot!(data2, xscale = :log10, yscale = :log10)
```

Now we want to fit a single powerlaw model to both of these spectra simultaneously, but with the powerlaw index fixed to be the same in both models.
```julia
model = PowerLaw()
prob = FittingProblem(model => data1, model => data2)

bind!(prob, :a)
prob
```

We can get a better look at our model configuration by using the [`details`](@ref) method:
```julia
details(prob)
```

In this printout we see that the `a` parameter of `Model 2` is bound to the `a` parameter of `Model 1`.

```julia
result = SpectralFitting.fit(prob, LevenbergMarquadt())

plot(data1, xscale = :log10, yscale = :log10)
plot!(data2, xscale = :log10, yscale = :log10)
plot!(result[1])
plot!(result[2])
```

Note that these fits are not perfect, because the underlying data have subtly different power law indices, but our fit is required to enforce the models to have the same value. If we release this requirement, the fit will be better, but the models will be entirely independent.

```julia
prob = FittingProblem(model => data1, model => data2)

result = SpectralFitting.fit(prob, LevenbergMarquadt())

plot(data1, xscale = :log10, yscale = :log10)
plot!(data2, xscale = :log10, yscale = :log10)
plot!(result[1])
plot!(result[2])
```