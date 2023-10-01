
@inline @generated function generated_model_call!(fluxes, energy, model, parameters)
    FunctionGeneration.generated_model_call!(fluxes, energy, model, parameters)
end

"""
    all_parameter_symbols(model::AbstractSpectralModel)

Returns a compile-time known tuple of all models symbols.
"""
@inline @generated function all_parameter_symbols(model::AbstractSpectralModel)
    params = FunctionGeneration.all_parameter_symbols(model)
    :($(params))
end

remake_with_parameters(model::AbstractSpectralModel, cache::ParameterCache) =
    _unsafe_remake_with_parameters(model, cache.parameters)
function remake_with_parameters(model::AbstractSpectralModel, params::AbstractArray)
    @assert length(params) == parameter_count(model)
    _unsafe_remake_with_parameters(model, params)
end
@inline @generated function _unsafe_remake_with_parameters(
    model::AbstractSpectralModel,
    parameters::AbstractArray,
)
    constructor =
        FunctionGeneration._construct_model_from_parameter_vector(model, parameters)
    :($(constructor))
end

"""
    objective_cache_count(model::AbstractSpectralModel)

Returns the number of flux arrays the model needs when using [`invokemodel!`](@ref).

# Example

```julia
model = XS_PhotoelectricAbsorption() * XS_PowerLaw()
objective_cache_count(model)
```
"""
@inline @generated function objective_cache_count(model::AbstractSpectralModel)
    FunctionGeneration.generated_maximum_flux_count(model)
end

@inline @generated function all_parameters_to_named_tuple(model::AbstractSpectralModel)
    FunctionGeneration.all_parameters_to_named_tuple(model)
end

@inline @generated function all_parameters_to_named_tuple(
    params::AbstractVector,
    model::AbstractSpectralModel,
)
    FunctionGeneration.all_parameters_to_named_tuple(params, model)
end

@inline @generated function parameter_count(model::AbstractSpectralModel)::Int
    params = FunctionGeneration.all_parameter_symbols(model)
    N = length(params)
    :($(N))
end

@inline @generated function parameter_count(model::CompositeModel)::Int
    info = SpectralFitting.FunctionGeneration.getinfo(model)
    N = reduce((total, i) -> (total + length(i.symbols)), info; init = 0)
    :($(N))
end

@inline @generated function model_parameters_tuple(model::AbstractSpectralModel)
    params = FunctionGeneration.model_parameters_tuple(model)
    :(($(params...),))
end

@inline @generated function _destructure_for_printing(model::CompositeModel)
    FunctionGeneration._destructure_for_printing(model)
end

export parameter_count
