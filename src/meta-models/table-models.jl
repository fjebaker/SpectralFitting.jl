"""
    abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{T,K} end

Abstract type representing table models, i.e. those models that interpolate or
load data from a table. 

First field in the struct **must be** `table`. See
[`PhotoelectricAbsorption`](@ref) for an example implementation.
"""
abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{T,K} end

# reflection tie-ins
Reflection.get_closure_symbols(::Type{<:AbstractTableModel}) = (:table,)

# `table` field is not a model parameter
Reflection.get_parameter_symbols(model::Type{<:AbstractTableModel}) =
    fieldnames(model)[2:end]

"""
    Base.copy(m::AbstractTableModel)

Create a copy of an [`AbstractTableModel`](@ref). This will copy all fields except
the `table` field, which is assumed to be a constant set of values that can be
shared by multiple copies.

When this is not the case, the user should redefine `Base.copy` for their particular
table model to copy the table as needed.
"""
function Base.copy(m::AbstractTableModel)
    typeof(m)(m.table, (copy(getproperty(m, f)) for f in fieldnames(typeof(m))[2:end])...)
end

abstract type AbstractCachedModel{T,K} <: AbstractSpectralModel{T,K} end

# reflection tie-ins
function Reflection.get_closure_symbols(::Type{<:AbstractCachedModel})
    (:cache,)
end
Reflection.get_parameter_symbols(model::Type{<:AbstractCachedModel}) =
    fieldnames(model)[2:end]

export AbstractTableModel
