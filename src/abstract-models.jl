export AbstractSpectralModel,
    AbstractSpectralModelKind,
    Multiplicative,
    Additive,
    Convolutional,
    modelkind,
    AbstractSpectralModelImplementation,
    XSPECImplementation,
    JuliaImplementation,
    implementation,
    AbstractSpectralModelClosureType,
    WithClosures,
    WithoutClosures,
    closurekind,
    has_closure_params,
    get_param,
    get_param_count,
    get_params,
    get_params_value,
    get_param_symbol_pairs,
    invokemodel,
    invokemodel!,
    flux_count,
    update_params!,
    freeze_parameter,
    free_parameter

"""
    abstract type AbstractSpectralModel{K}

Supertype of all spectral models. Sub-types must implement the following interface
- [`modelkind`](@ref)
- [`SpectralFitting.invoke!`](@ref)

The parametric type parameter `K` defines the [`AbstractSpectralModelKind`](@ref).
"""
abstract type AbstractSpectralModel{K} end
numbertype(::AbstractSpectralModel) = Sys.WORD_SIZE == 64 ? Float64 : Float32

"""
    abstract type AbstractSpectralModelKind

Abstract type of all model kinds. The algebra of models is as follows
```julia
A + A = A
M * M = M
M * A = A
C(A)  = A
```
where `A` is [`Additive`](@ref), `M` is [`Multiplicative`](@ref), and `C` is [`Convolutional`](@ref).
All other operations are prohibited, e.g. `C(M)` or `M * C`. To obtain `M * C` there must be an
additive component, e.g. `M * C(A)`.
"""
abstract type AbstractSpectralModelKind end
"""
    Additive <: AbstractSpectralModelKind
    Additive()

Additive models are effectively the sources of photons, and are the principle building blocks
of composite models. Every additive model has a normalisation parameter which re-scales the
flux by a constant factor `K`.

!!! note
    Defining custom additive models requires special care. See [Defining new models](@ref).
"""
struct Additive <: AbstractSpectralModelKind end
"""
    Multiplicative <: AbstractSpectralModelKind
    Multiplicative()

Multiplicative models act on [`Additive`](@ref) models, by element-wise
multiplying the flux in each energy bin of the additive model by a different factor.
"""
struct Multiplicative <: AbstractSpectralModelKind end
"""
    Convolutional <: AbstractSpectralModelKind
    Convolutional()

Convolutional models act on the flux generated by [`Additive`](@ref) models, similar to
[`Multiplicative`](@ref) models, however may convolve kernels through the flux also.
"""
struct Convolutional <: AbstractSpectralModelKind end

"""
    modelkind(M::Type{<:AbstractSpectralModel})
    modelkind(::AbstractSpectralModel)

Return the kind of model given by `M`: either `Additive`, `Multiplicative`, or `Convolutional`.
"""
modelkind(::Type{<:AbstractSpectralModel{K}}) where {K} = K()
modelkind(::M) where {M<:AbstractSpectralModel} = modelkind(M)

abstract type AbstractSpectralModelImplementation end
struct XSPECImplementation <: AbstractSpectralModelImplementation end
struct JuliaImplementation <: AbstractSpectralModelImplementation end

implementation(::Type{<:AbstractSpectralModel}) = JuliaImplementation()

abstract type AbstractSpectralModelClosureType end
struct WithClosures <: AbstractSpectralModelClosureType end
struct WithoutClosures <: AbstractSpectralModelClosureType end

closurekind(::Type{<:AbstractSpectralModel}) = WithoutClosures()

has_closure_params(::WithClosures) = true
has_closure_params(::WithoutClosures) = false
has_closure_params(M::Type{<:AbstractSpectralModel}) = has_closure_params(closurekind(M))
has_closure_params(::M) where {M<:AbstractSpectralModel} = has_closure_params(M)

# implementation interface
# never to be called directly
# favour `invokemodel!` instead
"""
    SpectralFitting.invoke!(flux, energy, M::Type{<:AbstractSpectralModel}, params...)

Used to define the behaviour of models. Should calculate flux of the model and write in-place
into `flux`.

!!! warning
    This function should not be called directly. Use [`invokemodel`](@ref) instead.

Parameters are passed in in-order as defined in the model structure. For example
```julia
struct MyModel{F1,F2,F3,...} <: AbstractSpectralModel
    p1::F1
    p2::F2
    p3::F3
    # ...
end
```
would have the arguments passed to `invoke!` as
```julia
function SpectralFitting.invoke!(flux, energy, ::Type{<:MyModel}, p1, p2, p3, ...)
    # ...
end
```

The only exception to this are [`Additive`](@ref) models, where the normalisation parameter
`K` is not passed to `invoke!`.
"""
invoke!(flux, energy, M::AbstractSpectralModel, params...) = error("Not defined for $(M).")

"""
    get_param(model::AbstractSpectralModel, s::Symbol)

Get a parameter from an [`AbstractSpectralModel`](@ref) by symbol.

# Example

```julia
model = XS_BlackBody()
get_param(model, :K)
```
"""
get_param(m::AbstractSpectralModel, s::Symbol) = getproperty(m, s)

"""
    get_param_count(model::AbstractSpectralModel)
    get_param_count(M::Type{<:AbstractSpectralModel})

Get the number of parameters a given [`AbstractSpectralModel`](@ref) has.

# Example

```julia
model = XS_BlackBody() + XS_PowerLaw()
get_param_count(model)
```
"""
get_param_count(::M) where {M<:AbstractSpectralModel} = get_param_count(M)
get_param_count(M::Type{<:AbstractSpectralModel}) = length(get_param_types(M))

"""
    get_params(m::AbstractSpectralModel)

Get a generator of all model parameters of an [`AbstractSpectralModel`](@ref).

# Example

```julia
model = XS_BlackBody() + XS_PowerLaw()
get_params(model)
```
"""
get_params(m::M) where {M<:AbstractSpectralModel} = (get_param(m, p) for p in FunctionGeneration.all_parameter_symbols(M))

"""
    get_params_value(m::AbstractSpectralModel)

Get a generator of all model parameter values of an [`AbstractSpectralModel`](@ref). See [`get_value`](@ref).

# Example

```julia
model = XS_BlackBody() + XS_PowerLaw()
get_params_value(model)
```
"""
get_params_value(m::AbstractSpectralModel) = (get_value(i) for i in get_params(m))
# todo: make this a proper iterator? also better name

"""
    get_param_symbol_pairs(m::AbstractSpectralModel)

Get a generator yielding `::Pair{Symbol,T}` of all model parameter values and their symbols,
for an [`AbstractSpectralModel`](@ref).

# Example

```julia
model = XS_BlackBody() + XS_PowerLaw()
get_param_symbol_pairs(model)
```
"""
get_param_symbol_pairs(m::M) where {M<:AbstractSpectralModel} =
    (p => get_param(m, p) for p in FunctionGeneration.all_parameter_symbols(M))

"""
    invokemodel(energy, model)
    invokemodel(energy, model, free_params)

Invoke the [`AbstractSpectralModel`](@ref) given by `model`, optionally overriding the free
parameters with values given in `free_params`. `free_params` may be a vector or tuple with element
type [`AbstractFitParameter`](@ref) or `Number`.

This function, unlike [`SpectralFitting.invoke!`](@ref) used to define models, is sensitive to performing
any normalisation or post-processing tasks that a specific model kind may require.

!!! note
    Users should always call models using [`invokemodel`](@ref) or [`invokemodel!`](@ref) to ensure
    normalisations and closures are accounted for.

`invokemodel` allocates the needed flux arrays based on the element type of `free_params` to allow
automatic differentation libraries to calculate parameter gradients.

In-place non-allocating variants are the [`invokemodel!`](@ref) functions.

# Example

```julia
model = XS_PowerLaw()
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, model)

p0 = [0.1, 2.0] # change K and a
invokemodel(energy, model, p0)
```
"""
function invokemodel(e, m::AbstractSpectralModel)
    flux = make_flux(eltype(e), length(e) - 1)
    invokemodel!(flux, e, m)
    flux
end
function invokemodel(e, m::AbstractSpectralModel, free_params)
    if eltype(free_params) <: Number
        # for compatability with AD
        flux = make_flux(eltype(free_params), length(e) - 1)
        invokemodel!(flux, e, m, free_params)
    else
        p0 = get_value.(free_params)
        flux = make_flux(eltype(p0), length(e) - 1)
        invokemodel!(flux, e, m, p0)
    end
    flux
end

"""
    invokemodel!(flux, energy, model)
    invokemodel!(flux, energy, model, free_params)
    invokemodel!(flux, energy, model, free_params, frozen_params)

In-place variant of [`invokemodel`](@ref), calculating the flux of an [`AbstractSpectralModel`](@ref)
given by `model`, optionally overriding the free and/or frozen parameter values. These arguments
may be a vector or tuple with element type [`AbstractFitParameter`](@ref) or `Number`.

The number of fluxes to allocate for a model may change if using any [`CompositeModel`](@ref)
as the `model`. It is generally recommended to use [`flux_count`](@ref) to ensure the correct number
of flux arrays are allocated with [`make_fluxes`](@ref) when using composite models.

Single spectral model components should use [`make_flux`](@ref) instead.

# Example

```julia
model = XS_PowerLaw()
energy = collect(range(0.1, 20.0, 100))
flux = make_flux(energy)
invokemodel!(flux, energy, model)

p0 = [0.1, 2.0] # change K and a
invokemodel!(flux, energy, model, p0)
```
"""
function invokemodel!(f, e, m::M) where {M<:AbstractSpectralModel}
    invokemodel!(f, e, M, get_params_value(m)...)
end
function invokemodel!(f, e, ::M, free_params) where {M<:AbstractSpectralModel}
    invokemodel!(f, e, M, free_params...)
end
# todo: do we want to support this function ?
# function invokemodel!(f, e, model::AbstractSpectralModel, free_params, frozen_params)
#     generated_model_call!(f, e, model, free_params, frozen_params)
# end
"""
    invokemodel!(flux, energy, ::Type{<:AbstractSpectralModel}, params...)

Specialisation of [`invokemodel!`](@ref) which dispatches on types. This is primarily used in
the generated functions, to assemble efficient model calls at compile time. In general, users
should prefer one of the other dispatches.
"""
invokemodel!(f, e, ::Type{M}, p...) where {M<:AbstractSpectralModel} =
    invokemodel!(f, e, modelkind(M), M, p...)
# implementations
@fastmath function invokemodel!(
    flux,
    energy,
    ::Additive,
    M::Type{<:AbstractSpectralModel},
    K,
    p...,
)
    invoke!(flux, energy, M, p...)
    flux .*= K
end
@fastmath function invokemodel!(
    flux,
    energy,
    ::AbstractSpectralModelKind,
    M::Type{<:AbstractSpectralModel},
    p...,
)
    invoke!(flux, energy, M, p...)
    flux
end

# bindings to generated functions

"""
    flux_count(model::AbstractSpectralModel)

Returns the number of flux arrays the model needs when using [`invokemodel!`](@ref).

# Example

```julia
model = XS_PhotoelectricAbsorption() * XS_PowerLaw()
flux_count(model)
```
"""
flux_count(model::AbstractSpectralModel) = generated_maximum_flux_count(model)


"""
    update_params!(model::AbstractSpectralModel, values)
    update_params!(model::AbstractSpectralModel, values, errors)

Update the free model parameters with `values`, optionally also updating the `errors`.
"""
function update_params!(model::AbstractSpectralModel, values)
    for (p, v, e) in zip(get_free_model_params(model), values)
        set_value!(p, v)
    end
    model
end
function update_params!(model::AbstractSpectralModel, values, errors)
    for (p, v, e) in zip(get_free_model_params(model), values, errors)
        set_value!(p, v)
        set_error!(p, e)
    end
    model
end


# printing

function modelinfo(m::M) where {M<:AbstractSpectralModel}
    params = join([get_value(p) for p in get_params(m)], ", ")
    "$(FunctionGeneration.model_base_name(M))[$(params)]"
end

function Base.show(io::IO, ::MIME"text/plain", m::M) where {M<:AbstractSpectralModel}
    params = [String(s) => p for (s, p) in get_param_symbol_pairs(m)]
    print(io, "$(FunctionGeneration.model_base_name(M))\n")

    pad = maximum(i -> length(first(i)), params) + 1

    for (s, val) in params
        print(io, "   $(rpad(s, pad)) => ")
        println(io, val)
    end
end

# parameter utilities
function remake_model(T::Type, model::AbstractSpectralModel)
    T(get_param_symbol_pairs(model)...)
end

function freeze_parameter(model, symbols...)
    error("Not implemented yet.")
end
function free_parameter(model, symbols...)
    error("Not implemented yet.")
end
