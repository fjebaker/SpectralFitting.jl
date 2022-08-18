# for future use: mission specific parsing
abstract type AbstractMission end

# abstract supertype
abstract type AbstractMetadata end

missiontrait(T::Type{<:AbstractMetadata}) = error("No trait set for $(T)")
missiontrait(::M) where {M<:AbstractMetadata} = missiontrait(M)

# interface
trim_meta!(::AbstractMetadata, inds) = nothing

group_meta(::AbstractMetadata, mask, inds) = nothing

# could be Response or Redistribution : how do we track this? 
struct ResponseMatrix{T}
    matrix::SparseMatrixCSC{T,Int}
    channels::Vector{Int}
    channel_energy_bins_low::Vector{T}
    channel_energy_bins_high::Vector{T}
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
end

abstract type AbstractDataset{T,M} end

# concrete type
mutable struct SpectralDataset{T,M} <: AbstractDataset{T,M}
    # store high and low seperately
    # incase discontinuous dataset
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
    counts::Vector{T}
    countserror::Vector{T}
    channels::Vector{Int}
    mask::BitVector
    meta::M
    response::ResponseMatrix{T}
end


missiontrait(::SpectralDataset{T,M}) where {T,M} = missiontrait(M)
