# Walkthrough

!!! warning
    This walk through has not been fleshed out with the relevant astrophysical content yet (for example, whether a fit is good, what the different parameters mean, etc.), and so assumes some familarity with spectral fitting in general.

    It is also not yet complete, nor a faithful illustration of everything SpectralFitting.jl can do. It serves to illustrate similarities and differences in syntax between SpectralFitting.jl and XSPEC.

This example walkthrough is the SpectralFitting.jl equivalent of the [Walk through XSPEC](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/node35.html) from the XSPEC manual. We will use the same dataset, available for download from this [link to the data files](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/walkthrough.tar.gz).

## Overview

The first thing we want to do is load our datasets. Unlike in XSPEC, we have no requirement of being in the same directory as the data, or even that all of the response, ancillary, and spectral files be in the same place. For simplicity, we'll assume they are:

!!! note 
    Be sure to set `DATADIR` pointing to the directory where you keep the walkthrough data.

```@example walk
using SpectralFitting, Plots

DATADIR = "..."
DATADIR = length(get(ENV, "CI", "")) > 0 ? @__DIR__() * "/../../ex-datadir" : "/home/lilith/Developer/jl/datasets/xspec/walkthrough" # hide
spec1_path = joinpath(DATADIR, "s54405.pha")
data = OGIPDataset(spec1_path) 
```

This will print a little card about our data, which shows us what else SpectralFitting.jl loaded. We can see the `Primary Spectrum`, the `Response`, but that the `Background` and `Ancillary` response files are missing. That's to be expected, since we don't have those files in the dataset. 

We can check what paths it used by looking at
```@example walk
data.paths
```

We can load and alter any part of a dataset as we do our fitting. For example, if you have multiple different ancillary files at hand, switching them between fits is a one-liner.

To visualize our data, we can use some of the [Plots.jl](https://docs.juliaplots.org/latest/) recipes included in SpectralFitting.jl:

```@example walk
plot(data, xlims = (0.5, 70), xscale = :log10)
```

Note that the units are currently not divided by the energy bin widths. We can either do that manually, or use the [`normalize!`](@ref) to convert whatever units the data is currently in to the defacto standard `counts s⁻¹ keV⁻¹` for fitting. Whilst we're at it, we see in the model card that there are 40 [bad quality bins](https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007.pdf) still present in our data. We can drop those as well, and plot the data on log-log axes:

```@example walk
normalize!(data)
drop_bad_channels!(data)
plot(data, ylims = (0.001, 2.0), yscale = :log10, xscale = :log10)
```

Note that when there are no negative axes, the scale defaults to log on the plot unless otherwise specified.

Next we want to specify a model to fit to this data. Models that are prefixed with `XS_` are models that are linked from the XSPEC model library, provided via [LibXSPEC_jll](https://github.com/astro-group-bristol/LibXSPEC_jll.jl). For a full list of the models, see [Models library](@ref).

!!! warning
    It is advised to **use the Julia implemented models**. This allows various calculations to benefit from automatic differentiation, efficient multi-threading, GPU offloading, and various other useful things, see [Why & how](@ref).

We will start by fitting a photoelectric absorption model that acts on a power law model:

!!! note
    To see information about a model, use the `?` in the Julia REPL:
    ```julia
    julia> ?PowerLaw
    XS_PowerLaw(K, a)

        •  K: Normalisation.

        •  a: Photon index.

    Example
    ≡≡≡≡≡≡≡
    ...
    ```

```@example walk
model = PhotoelectricAbsorption() * PowerLaw()
```

If we want to specify paramters of our model at instantiation, we can do that with
```@example walk
model = PhotoelectricAbsorption() * PowerLaw(a = FitParam(3.0))
```

SpectralFitting.jl adopts the SciML problem-solver abstraction, so to fit a model to data we specify a [`FittingProblem`](@ref):

```@example walk
prob = FittingProblem(model => data)
```

SpectralFitting.jl makes a huge wealth of optimizers availble from [Optimizations.jl](https://github.com/SciML/Optimization.jl), and others from further afield. For consistency with XSPEC, we'll use here a delayed-gratification least-squares algorithm from [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl):

```@example walk
result = fit(prob, LevenbergMarquadt())
```

Here we can see the parameter vector, the estimated error on each parameter, and the measure of the fit statistic (here chi squared). We can overplot our result on our data easily:

```@example walk
plot(data, 
    ylims = (0.001, 2.0), 
    xscale = :log10, 
    yscale = :log10
)
plot!(data, result)
```

Our model does not account for the high energy range well. We can ignore that range for now, and select everything from 0 to 15 keV and refit:

```@example walk
mask_energies!(data, 0, 15)
result = fit(prob, LevenbergMarquadt())
```

```@example walk
plot(data, 
    ylims = (0.001, 2.0), 
    xscale = :log10, 
    yscale = :log10
)
plot!(data, result, label = "PowerLaw")
```

The result is not yet baked into our model, and represents just the outcome of the fit. To update the parameters and errors in the model, we can use [`update_model!`](@ref)

```@example walk
update_model!(model, result)
```

!!! note
    Since fitting and updating a model is often done in tandem, SpectralFitting.jl has both a [`fit`](@ref) and [`fit!`](@ref) method, the latter automatically updates the model parameters after fit.

To estimate the goodness of our fit, we can mimic the `goodness` command from XSPEC. This will use the [`simulate`](@ref) function to simulate spectra for a dataset (here determined by the result), and fit the model to the simulated dataset. The fit statistic for each fit is then appended to an array, which we can use to plot a histogram:

```@example walk
spread = goodness(result; N = 1000, seed = 42)
histogram(spread, ylims = (0, 300), label = "Simulated")
vline!([result.χ2], label = "Best fit")
```

Note we have set the random number generator seed with `seed = 42` to allow our results to be strictly reproduced.

The `goodness` command will log the percent of simulations with a fit statistic better than the result, but we can equivalently calculate that ourselves:

```@example walk
count(<(result.χ2), spread) * 100 / length(spread)
```

Next we want to calculate the flux in an energy range observed by the detector. We can do this with [`LogFlux`](@ref) or [`XS_CalculateFlux`](@ref), as they are both equivalent implementations.

We can modify our model by accessing properties from the model card and writing a new expression:

```@example walk
calc_flux = XS_CalculateFlux(
    E_min = FitParam(0.2, frozen = true), 
    E_max = FitParam(2.0, frozen = true),
    log10Flux = FitParam(-10.3, lower_limit = -100, upper_limit = 100),
)

flux_model = model.m1 * calc_flux(model.a1)
```

Since we used the old model to define the new one, our best fit values are automatically copied into the new model. We can now freeze the normalization, as we are using the flux integrating model to scale the powerlaw component:

```@example walk
flux_model.a1.K.frozen = true
flux_model
```

Looking at the data card, we see the fit domain does not include the full region that we want to integrate the flux over. We therefore need to extend the fitting domain:
```@example walk
flux_problem = FittingProblem(flux_model => data)
# TODO: domain extensions not fully implemented yet
```

Now to fit we can repeat the above procedure, and even overplot the region of flux we integrated:
```@example walk
flux_result = fit(flux_problem, LevenbergMarquadt())

plot(data, 
    ylims = (0.001, 2.0), 
    xscale = :log10, 
    yscale = :log10
)
plot!(data, flux_result)
vspan!([flux_model.c1.E_min.value, flux_model.c1.E_max.value], alpha = 0.5)
```

Let's try alternative models to see how they fit the data. First, an absorbed black body:

```@example walk
model2 = PhotoelectricAbsorption() * XS_BlackBody()
```

We fit in the same way as before:

```@example walk
prob2 = FittingProblem(model2 => data)
result2 = fit!(prob2, LevenbergMarquadt())
```

Let's overplot this result against our power law result:

```@example walk
dp = plot(data, 
    ylims = (0.001, 2.0), 
    xscale = :log10, 
    yscale = :log10,
    legend = :bottomleft,
)
plot!(dp, data, result, label = "PowerLaw $(round(result.χ2))")
plot!(dp, data, result2, label = "BlackBody $(round(result2.χ2))")
```

Or a bremsstrahlung model:

```@example walk
model3 = PhotoelectricAbsorption() * XS_BremsStrahlung()
prob3 = FittingProblem(model3 => data)
result3 = fit(prob3, LevenbergMarquadt())
```

```@example walk
plot!(dp, data, result3, label = "Brems $(round(result3.χ2))")
```

Let's take a look at the residuals of these three models. There are utility methods for this in SpectralFitting.jl, but we can easily just interact with the result directly:

```@example walk
function residuals(result)
    # select which result we want (only have one, but for generalisation to multi-model fits)
    r = result[1] 
    y = invoke_result(r)
    @. (r.objective - y) / sqrt(r.variance)
end

domain = SpectralFitting.plotting_domain(data)

rp = hline([0], linestyle = :dash, legend = false)
plot!(rp,domain, residuals(result), seriestype = :stepmid)
plot!(rp, domain, residuals(result2), seriestype = :stepmid)
plot!(rp, domain, residuals(result3), seriestype = :stepmid)
rp
```

We can compose this figure with our previous one, and change to a linear x scale:

```@example walk
plot(dp, rp, layout = grid(2, 1, heights = [0.7, 0.3]), link = :x, xscale = :linear)
```

Let's modify the black body model with a continuum component

```@example walk
bbpl_model = model2.m1 * (PowerLaw() + model2.a1) |> deepcopy
```

!!! note
    We pipe the model to `deepcopy` to create a copy of all the model parameters. Not doing this means the parameters in `bbpl_model` will be aliased to the parameters in `model2`, and changing one with change the other.

We'll freeze the hydrogen column density parameter to the galactic value and refit:

```@example walk
bbpl_model.ηH_1.value = 4
bbpl_model.ηH_1.frozen = true
bbpl_model
```

And fitting:

```@example walk
bbpl_result = fit(
    FittingProblem(bbpl_model => data), 
    LevenbergMarquadt()
)
```

Let's plot the result:

```@example walk
plot(data, 
    ylims = (0.001, 2.0), 
    xscale = :log10, 
    yscale = :log10,
    legend = :bottomleft,
)
plot!(data, bbpl_result)
```

Update the model and fix the black body temperature to 2 keV:

```@example walk
update_model!(bbpl_model, bbpl_result)

bbpl_model.T_1.value = 2.0
bbpl_model.T_1.frozen = true
bbpl_model
```

Fitting:

```@example walk
bbpl_result2 = fit(
    FittingProblem(bbpl_model => data), 
    LevenbergMarquadt()
)
```

Overplotting this new result:

```@example walk
plot!(data, bbpl_result2)
```