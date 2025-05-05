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
Base.copy(m::AsConvolution) =
    AsConvolution(copy(backing_model(m)), m.domain, deepcopy(m.cache))

function parameter_count(m::AsConvolution{<:AbstractSpectralModel{T,K}}) where {T,K}
    count = parameter_count(backing_model(m))
    if K <: Additive
        count - 1
    else
        count
    end
end

function remake_with_number_type(model::AsConvolution)
    AsConvolution(remake_with_number_type(backing_model(model)), model.domain, model.cache)
end

function remake_with_parameters(
    model::AsConvolution{<:AbstractSpectralModel{T,K}},
    params::Tuple,
) where {T,K}
    _params = if K <: Additive
        # Need an additional parameter for the normalisation term
        (params[1], params...)
    else
        params
    end
    AsConvolution(remake_with_parameters(model.model, _params), model.domain, model.cache)
end

_model_name(model::AsConvolution) = "AsConvolution[$(_model_name(model.model))]"

function unpack_parameters_as_named_tuple(
    model::AsConvolution{<:AbstractSpectralModel{T,K}},
) where {T,K}
    ps = unpack_parameters_as_named_tuple(model.model)
    if K <: Additive
        # TODO: drop the K
        NamedTuple{keys(ps)[2:end]}((ps...,)[2:end])
    else
        ps
    end
end

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
    invoke!(model_output, model.domain, backing_model(model))

    # do the convolution
    convolve!(convolution_cache, output, domain, model_output, model.domain)

    # overwrite the output
    @. output = convolution_cache
end

export AsConvolution
