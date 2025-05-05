mutable struct CacheEntry{T}
    cache::Vector{T}
    params::Vector{T}
    domain_limits::Tuple{T,T}
    size_of_element::Int

    function CacheEntry(params::AbstractVector{<:Number})
        T = eltype(params)
        cache = zeros(T, 1)
        new{T}(cache, params, (zero(T), zero(T)), sizeof(T))
    end
end

"""
    AutoCache

Used to automatically create a cache of another model, to avoid re-evaluating
the model if the next parameters are close to the previous parameters. The
intended use is for fitting expensive models which.

## Example

```julia
model = PhotoelectricAbsorption() * AutoCache(PowerLaw())
```

In the above model, the [`PowerLaw`](@ref) component will be augmented with the
caching behaviour.
"""
struct AutoCache{M,T,K,C<:CacheEntry} <: AbstractModelWrapper{M,T,K}
    model::M
    cache::C
    abstol::Float64
    enabled::Bool
    function AutoCache(
        model::AbstractSpectralModel{T,K},
        cache::CacheEntry,
        abstol,
        enabled::Bool,
    ) where {T,K}
        new{typeof(model),T,K,typeof(cache)}(model, cache, abstol, enabled)
    end
end

function Base.copy(m::AutoCache)
    AutoCache(copy(m.model), deepcopy(m.cache), m.abstol, m.enabled)
end

_model_name(model::AutoCache) = "AutoCache[$(_model_name(model.model))]"
unpack_parameters_as_named_tuple(model::AutoCache) =
    unpack_parameters_as_named_tuple(model.model)
remake_with_number_type(model::AutoCache) = AutoCache(
    remake_with_number_type(model.model),
    model.cache,
    model.abstol,
    model.enabled,
)

function AutoCache(
    model::AbstractSpectralModel{T,K};
    abstol = 1e-3,
    enabled = true,
) where {T,K}
    params = [get_value(p) for p in unpack_parameters_as_named_tuple(model)]
    cache = CacheEntry(params)
    AutoCache(model, cache, abstol, enabled)
end

function _reinterpret_dual(
    M::Type{<:AbstractSpectralModel},
    ::Type,
    v::AbstractArray,
    n::Int,
)
    needs_resize = n > length(v)
    if needs_resize
        @warn "$(Base.typename(M).name): Growing dual buffer..."
        resize!(v, n)
    end
    view(v, 1:n), needs_resize
end
function _reinterpret_dual(
    M::Type{<:AbstractSpectralModel},
    DualType::Type{<:ForwardDiff.Dual},
    v::AbstractArray{T},
    n::Int,
) where {T}
    n_elems = div(sizeof(DualType), sizeof(T)) * n
    needs_resize = n_elems > length(v)
    if needs_resize
        @warn "$(Base.typename(M).name): Growing dual buffer..."
        resize!(v, n_elems)
    end
    reinterpret(DualType, view(v, 1:n_elems)), needs_resize
end

function invoke!(output, domain, model::AutoCache{M,T,K}) where {M,T,K}
    if model.enabled == false
        return invoke!(output, domain, model.model)
    end
    D = promote_type(eltype(domain), T)

    _new_params = unpack_parameters_as_named_tuple(model.model)
    _new_limits = (first(domain), last(domain))

    output_cache, out_resized =
        _reinterpret_dual(typeof(model), D, model.cache.cache, length(output))
    param_cache, _ =
        _reinterpret_dual(typeof(model), D, model.cache.params, length(_new_params))

    same_domain = model.cache.domain_limits == _new_limits

    # if the parameter size has changed, need to rerun the model
    if (!out_resized) && (model.cache.size_of_element == sizeof(D)) && (same_domain)
        # ignore the normalisation, since that's applied later
        _pc, _np = @views if K === Additive
            param_cache[2:end], _new_params[2:end]
        else
            param_cache, _new_params
        end
        # if all parameters within some tolerance, then just return the cache
        within_tolerance = all(zip(_pc, _np)) do I
            p, pm = I
            abs((get_value(p) - get_value(pm)) / p) < model.abstol
        end

        if within_tolerance
            @. output = output_cache
            return output
        end
    end

    model.cache.size_of_element = sizeof(D)
    invoke!(output_cache, domain, model.model)
    # update the auto cache infos
    model.cache.domain_limits = _new_limits

    @. param_cache = get_value(_new_params)
    @. output = output_cache
end

export AutoCache
