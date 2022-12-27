abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{K} end

closurekind(::Type{<:AbstractTableModel}) = WithClosures()

get_closure_param_fields(::Type{<:AbstractTableModel}) = (:table,)
# drop the table parameter
get_param_types(T::Type{<:AbstractTableModel}) = T.types[2:end]
get_param_symbols(T::Type{<:AbstractTableModel}) = fieldnames(T)[2:end]

# runtime access
get_param_symbols(m::T) where {T<:AbstractTableModel} = get_param_symbols(T)

function invokemodel!(f, e, m::M) where {M<:AbstractTableModel}
    invokemodel!(f, e, M, m.table, get_params_value(m)...)
end

# todo: use subtypes to ensure everything is correct in the model definitions
