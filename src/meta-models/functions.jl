"""
    AsConvolution

Turn an additive model into a convolutional model.

## Example

```
convolution_model = AsConvolution(GaussianLine())
```

The above model will now convolve the [`GaussianLine`](@ref) model onto whatever
it is applied to.
"""
struct AsConvolution{M,T,V,P} <: AbstractModelWrapper{M,T,Convolutional}
    model::M
    # the domain on which we evaluate this model
    domain::V
    # an additional output cache
    cache::NTuple{2,Vector{P}}
    function AsConvolution(
        model::AbstractSpectralModel{T},
        domain::V,
        cache::NTuple{2,Vector{P}},
    ) where {T,V,P}
        new{typeof(model),T,V,P}(model, domain, cache)
    end
end

"""
    Base.copy(m::AsConvolution)

Creates a copy of an [`AsConvolution`](@ref) wrapped model. Will make a 
`deepcopy` of the cache to elimiate possible thread contention, but does not
copy the domain.
"""
Base.copy(m::AsConvolution) = AsConvolution(copy(m.model), m.domain, deepcopy(m.cache))

function AsConvolution(
    model::AbstractSpectralModel{T};
    domain = collect(range(0, 2, 100)),
) where {T}
    output = invokemodel(domain, model)
    AsConvolution(model, domain, (output, deepcopy(output)))
end

function invoke!(output, domain, model::AsConvolution{M,T}) where {M,T}
    D = promote_type(eltype(domain), T)
    model_output, _ =
        _reinterpret_dual(typeof(model), D, model.cache[1], length(model.domain) - 1)
    convolution_cache, _ =
        _reinterpret_dual(typeof(model), D, model.cache[2], length(output))

    # invoke the child model
    invoke!(model_output, model.domain, model.model)

    # do the convolution
    convolve!(convolution_cache, output, domain, model_output, model.domain)

    # overwrite the output
    @views output .= convolution_cache
end

function Reflection.get_parameter_symbols(
    ::Type{<:AsConvolution{M}},
) where {M<:AbstractSpectralModel{T,K}} where {T,K}
    syms = Reflection.get_parameter_symbols(M)
    if K === Additive
        # we need to lose the normalisation parameter
        (syms[2:end]...,)
    else
        syms
    end
end

function Reflection.make_constructor(
    M::Type{<:AsConvolution{Model}},
    closures::Vector,
    params::Vector,
    T::Type,
) where {Model<:AbstractSpectralModel{Q,K}} where {Q,K}
    num_closures = fieldcount(M) - 1 # ignore the `model` field
    my_closures = closures[1:num_closures]

    model_params = if K === Additive
        # insert a dummy normalisation to the constructor
        vcat(:(one($T)), params)
    else
        params
    end

    model_constructor =
        Reflection.make_constructor(Model, closures[(num_closures+1):end], model_params, T)
    :($(Base.typename(M).name)($(model_constructor), $(my_closures...)))
end

export AsConvolution