# these should work for both single models and composite models
# so that the parameters can all be allocated in one go
# for multi model: give a view to each cache
struct ParameterCache{M<:AbstractArray,V,T<:Number}
    free_mask::M # bit vector or a view into one
    parameters::V
    frozen_values::Vector{T}
end

function _make_free_mask(params::AbstractArray{<:FitParam})
    free_mask = BitVector(undef, length(params))
    for (i, p) in enumerate(params)
        free_mask[i] = isfree(p)
    end
    free_mask
end

function ParameterCache(params::AbstractArray{<:FitParam})
    free_mask = _make_free_mask(params)
    frozen = params[.!free_mask]
    ParameterCache(free_mask, map(get_value, params), map(get_value, frozen))
end

function _update_conditional!(parameters, mask, new_parameters, frozen)
    j::Int = 1
    k::Int = 1
    for (i, free) in enumerate(mask)
        if free
            parameters[i] = new_parameters[j]
            j += 1
        else
            parameters[i] = frozen[k]
            k += 1
        end
    end
end

_get_parameters(cache::ParameterCache, params) = cache.parameters
_get_parameters(cache::ParameterCache{M,V}, params) where {M<:AbstractArray,V<:DiffCache} =
    get_tmp(cache.parameters, params)

function update_free_parameters!(cache::ParameterCache, params)
    @assert count(cache.free_mask) == length(params)
    _update_conditional!(
        _get_parameters(cache, params),
        cache.free_mask,
        params,
        cache.frozen_values,
    )
    cache
end

export ParameterCache, update_free_parameters!
