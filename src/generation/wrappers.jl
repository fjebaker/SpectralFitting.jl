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


"""
    all_parameter_symbols(model::AbstractSpectralModel)

Returns a compile-time known tuple of all models symbols.
This method is not defined for [`CompositeModel`](@ref). Prefer [`modelparameters`](@ref).
"""
@generated function all_parameter_symbols(model::AbstractSpectralModel)
    params = FunctionGeneration.all_parameter_symbols(model)
    :($(params))
end
all_parameter_symbols(::CompositeModel) = throw("This inspection method is for base models only.")


"""
    free_parameter_symbols(model::AbstractSpectralModel)

Returns a compile-time known tuple of symbols corresponding to those parameters which are free in the model.
This method is not defined for [`CompositeModel`](@ref). Prefer [`modelparameters`](@ref).
"""
@generated function free_parameter_symbols(model::AbstractSpectralModel)
    params = FunctionGeneration.free_parameter_symbols(model)
    :($(params))
end
free_parameter_symbols(::CompositeModel) = throw("This inspection method is for base models only.")

"""
    frozen_parameter_symbols(model::AbstractSpectralModel)

Returns a compile-time known tuple of symbols corresponding to those parameters which are frozen in the model.
This method is not defined for [`CompositeModel`](@ref). Prefer [`modelparameters`](@ref).
"""
@generated function frozen_parameter_symbols(model::AbstractSpectralModel)
    params = FunctionGeneration.frozen_parameter_symbols(model)
    :($(params))
end
frozen_parameter_symbols(::CompositeModel) = throw("This inspection method is for base models only.")
