# these should work for both single models and composite models
# so that the parameters can all be allocated in one go
# for multi model: give a view to each cache
struct ParameterCache{M<:AbstractArray,V}
    free_mask::M # bit vector or a view into one
    parameters::V
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
    ParameterCache(free_mask, map(get_value, params))
end

function _update_conditional!(parameters, mask, new_parameters, condition)
    j::Int = 1
    for (i, free) in enumerate(mask)
        if condition(free)
            parameters[i] = new_parameters[j]
            j += 1
        end
    end
end

_get_parameters(cache::ParameterCache, params) = cache.parameters
_get_parameters(cache::ParameterCache{M,V}, params) where {M<:AbstractArray,V<:DiffCache} =
    get_tmp(cache.parameters, params)

function update_free_parameters!(cache::ParameterCache, params)
    @assert count(cache.free_mask) == length(params)
    _update_conditional!(_get_parameters(cache, params), cache.free_mask, params, ==(true))
    cache
end

function update_frozen_parameters!(cache::ParameterCache, params)
    parameters = _get_parameters(cache, params)
    @assert length(parameters) - count(cache.free_mask) == length(params)
    _update_conditional!(parameters, cache.free_mask, params, ==(false))
    cache
end

export ParameterCache, update_free_parameters!, update_frozen_parameters!
