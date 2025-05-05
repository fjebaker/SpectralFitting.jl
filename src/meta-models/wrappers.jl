"""
    abstract type AbstractModelWrapper{M,T,K} <: AbstractSpectralModel{T,K} end

First field of the struct must be `model`.
"""
abstract type AbstractModelWrapper{M<:AbstractSpectralModel,T,K} <:
              AbstractSpectralModel{T,K} end

# tie in dispatches              
backing_model(model::AbstractModelWrapper) = getfield(model, :model)
closure_and_parameter(model::AbstractModelWrapper) =
    closure_and_parameter(backing_model(model))

normalisation(model::AbstractModelWrapper{M,T,Additive}) where {M,T} =
    normalisation(backing_model(model))

function remake_with_number_type(m::AbstractModelWrapper)
    params = (unpack_parameters_as_named_tuple(backing_model(m))...,)
    remake_with_parameters(m, ((get_value(p) for p in params)...,))
end

function Base.propertynames(m::AbstractModelWrapper)
    keys(unpack_parameters_as_named_tuple(backing_model(m)))
end

function Base.getproperty(m::AbstractModelWrapper{<:Any,<:FitParam}, symb::Symbol)
    params = unpack_parameters_as_named_tuple(backing_model(m))
    if haskey(params, symb)
        params[symb]
    else
        getfield(m, symb)
    end
end
