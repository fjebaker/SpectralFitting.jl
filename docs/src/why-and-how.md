# Why & how

SpectralFitting.jl is a package for fitting models to spectral data, similar to [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/), [Sherpa](https://sherpa.readthedocs.io/en/latest/ciao.html) or [ISIS](https://space.mit.edu/CXC/isis/).

The rationale for this package is to provide a unanimous interface for different model libraries, and to leverage the bleeding edge advancements in computing that are available in Julia, including the rich statistics ecosystem, with [automatic-differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) and [_speed_](https://julialang.org/benchmarks/).

SpectralFitting aims to provide highly optimised and flexible fitting algorithms, along with a library of spectral models, for use in any field of Astronomy that concerns itself with spectral data.

## Rewriting model calls during invocation

```@setup model_invocation
using SpectralFitting
```

SpectralFitting.jl tries to optimise model invocation through source-rewriting. For compatibility with the XSPEC model library, this is achieved through aggressive pre-allocation and shuffling of output vectors. All XSPEC models output the result of their calculation through side effects into a flux array passed as an argument, and therefore each model invocation requires its own output before addition or multiplication of fluxes may occur.

Principally, combining several models together would look like this:

```@example model_invocation
energy = collect(range(0.1, 20.0, 100))

flux1 = invokemodel(energy, XS_PowerLaw())
flux2 = invokemodel(energy, XS_PowerLaw(a=FitParam(3.0)))
flux3 = invokemodel(energy, XS_BlackBody())
flux4 = invokemodel(energy, XS_PhotoelectricAbsorption())

total_flux = @. flux4 * (flux1 + flux2 + flux3)
sum(total_flux)
```

But these operations could also be performed in a different order:

```@example model_invocation
flux1 = invokemodel(energy, XS_PowerLaw())
flux2 = invokemodel(energy, XS_PowerLaw(a=FitParam(3.0)))
total_flux = @. flux1 + flux2

flux3 = invokemodel(energy, XS_BlackBody())
@. total_flux = total_flux + flux3

flux4 = invokemodel(energy, XS_PhotoelectricAbsorption())
@. total_flux = total_flux * flux4
sum(total_flux)
```

Doing so would allow us to only pre-allocate 2 flux arrays, instead of 4 when using the in-place variants:

```@example model_invocation
flux1, flux2 = make_fluxes(energy, 2)

invokemodel!(flux1, energy, XS_PowerLaw())
invokemodel!(flux2, energy, XS_PowerLaw(a=FitParam(3.0)))
@. flux1 = flux1 + flux2

invokemodel!(flux2, energy, XS_BlackBody())
@. flux1 = flux1 + flux2

invokemodel!(flux2, energy, XS_PhotoelectricAbsorption())
@. flux1 = flux1 * flux2
sum(flux1)
```

It is precisely this re-writing that SpectralFitting performs via [`@generated`](https://docs.julialang.org/en/v1/manual/metaprogramming/#Generated-functions) functions. We can inspect the code used to generate the invocation body after defining a [`CompositeModel`](@ref):

```@example model_invocation
fluxes = (flux1, flux2)

model = XS_PhotoelectricAbsorption() * (
    XS_PowerLaw() + XS_PowerLaw(a=FitParam(3.0)) + XS_BlackBody()
)

params = get_value.(get_params(model))

SpectralFitting.__generated_model_call!(fluxes, energy, typeof(model), params)
nothing # hide
```

```@raw html
<pre class="documenter-example-output"><code class="nohighlight hljs ansi">begin
    @inbounds let (flux1, flux2) = fluxes
        var"##K#356" = params[1]
        var"##a#357" = params[2]
        var"##K#358" = params[3]
        var"##a#359" = params[4]
        var"##K#360" = params[5]
        var"##T#361" = params[6]
        var"##ηH#362" = params[7]
        invokemodel!(flux1, energy, XS_PowerLaw, var"##K#356", var"##a#357")
        invokemodel!(flux2, energy, XS_PowerLaw, var"##K#358", var"##a#359")
        @. flux1 = flux1 + flux2
        invokemodel!(flux2, energy, XS_BlackBody, var"##K#360", var"##T#361")
        @. flux1 = flux1 + flux2
        invokemodel!(flux2, energy, XS_PhotoelectricAbsorption, var"##ηH#362")
        @. flux1 = flux1 * flux2
        return flux1
    end
end
</code><button class="copy-button fas fa-copy"></button></pre>
```

This generated function also takes care of some other things for us, such as unpacking parameters (optionally unpacking frozen parameters separately), and ensuring any closure are passed to [`invokemodel`](@ref) if a model needs them (e.g., [`SurrogateSpectralModel`](@ref)).

This is achieved by moving as much information as possible about the model and its construction to its type, such that all of the invocation and parameter unpacking may be inferred at compile time.

Naturally, the [`CompositeModel`](@ref) types also support the out-of-place [`invokemodel`](@ref) and will allocate the minimum number of flux arrays needed, inferred using [`flux_count`](@ref):

```@example model_invocation
flux = invokemodel(energy, model)
sum(flux)
```

!!! note
    With the addition of more pure-Julia models, non-allocating methods without aggressive pre-allocation are possible, and will be added in the future. Such methods may allow models to add or multiply in-place on the total flux array, instead of relying on later broadcasts.