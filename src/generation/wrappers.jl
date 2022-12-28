@generated function generated_model_parameter_count(model)
    :($(FunctionGeneration.model_parameter_count(model)))
end

@generated function generated_frozen_model_parameter_count(model)
    :($(FunctionGeneration.model_frozen_parameter_count(model)))
end

@generated function generated_free_model_parameter_count(model)
    :($(FunctionGeneration.model_free_parameter_count(model)))
end

@generated function generated_model_call!(fluxes, energy, model, params)
    FunctionGeneration.generated_model_call!(fluxes, energy, model, params)
end

@generated function generated_model_call!(fluxes, energy, model, free_params, frozen_params)
    FunctionGeneration.generated_model_call!(
        fluxes,
        energy,
        model,
        free_params,
        frozen_params,
    )
end

@generated function generated_maximum_flux_count(model)
    FunctionGeneration.generated_maximum_flux_count(model)
end

@generated function generated_model_types(model)
    ga = assemble_aggregate_info(model)
    :(($(ga.models...),))
end

@generated function generated_model_parameter_type(model)
    T = FunctionGeneration.model_T(model)
    :($(T))
end

"""
    all_parameter_symbols(model::AbstractSpectralModel)

Returns a compile-time known tuple of all models symbols.
"""
@generated function all_parameter_symbols(model::AbstractSpectralModel)
    params = FunctionGeneration.all_parameter_symbols(model)
    :($params)
end

"""
    free_parameter_symbols(model::AbstractSpectralModel)

Returns a compile-time known tuple of symbols corresponding to those parameters which are free in the model.
"""
@generated function free_parameter_symbols(model::AbstractSpectralModel)
    params = FunctionGeneration.free_parameter_symbols(model)
    :($params)
end

"""
    frozen_parameter_symbols(model::AbstractSpectralModel)

Returns a compile-time known tuple of symbols corresponding to those parameters which are frozen in the model.
"""
@generated function frozen_parameter_symbols(model::AbstractSpectralModel)
    params = FunctionGeneration.frozen_parameter_symbols(model)
    :($params)
end
