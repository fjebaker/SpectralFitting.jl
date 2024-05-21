"""
    abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{T,K} end

Abstract type representing table models, i.e. those models that interpolate or
load data from a table. 

First field in the struct **must be** `table`. See
[`PhotoelectricAbsorption`](@ref) for an example implementation.
"""
abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{T,K} end

# reflection tie-ins
function Reflection.get_closure_symbols(::Type{<:AbstractTableModel})
    (:table,)
end
function Reflection.get_parameter_symbols(model::Type{<:AbstractTableModel})
    # `table` field is not a model parameter
    fieldnames(model)[2:end]
end

export AbstractTableModel
