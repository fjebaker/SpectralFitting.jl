abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{K} end

closurekind(::Type{<:AbstractTableModel}) = WithClosures()

function FunctionGeneration.closure_parameter_symbols(::Type{<:AbstractTableModel})
    (:table,)
end
function FunctionGeneration.all_parameter_symbols(model::Type{<:AbstractTableModel})
    # `table` field is not a model parameter
    fieldnames(model)[2:end]
end

function invokemodel!(f, e, m::AbstractTableModel{T,K}) where {T,K}
    # pass table as the last argument
    invokemodel!(f, e, typeof(m), get_params_value(m)..., m.table)
end

# todo: use subtypes to ensure everything is correct in the model definitions

export AbstractTableModel
