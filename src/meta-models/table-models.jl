abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{T,K} end

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
    M{typeof(model.table),T}(model.table, new_params...)
end

function Reflection.get_closure_symbols(::Type{<:AbstractTableModel})
    (:table,)
end
function Reflection.get_parameter_symbols(model::Type{<:AbstractTableModel})
    # `table` field is not a model parameter
    fieldnames(model)[2:end]
end

export AbstractTableModel
