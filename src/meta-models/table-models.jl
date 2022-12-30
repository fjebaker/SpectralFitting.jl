abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{T,K} end

closurekind(::Type{<:AbstractTableModel}) = WithClosures()

function FunctionGeneration.closure_parameter_symbols(::Type{<:AbstractTableModel})
    (:table,)
end
function FunctionGeneration.all_parameter_symbols(model::Type{<:AbstractTableModel})
    # `table` field is not a model parameter
    fieldnames(model)[2:end]
end

function remake_with_number_type(model::AbstractTableModel{FitParam{T}}) where {T}
    M = typeof(model).name.wrapper
    params = modelparameters(model)
    new_params = convert.(T, params)
    M{typeof(model.table),T,modelkind(typeof(model))}(model.table, new_params...)
end

# todo: use subtypes to ensure everything is correct in the model definitions

export AbstractTableModel
