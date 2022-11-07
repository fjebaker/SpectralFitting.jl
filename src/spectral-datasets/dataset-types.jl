# for future use: mission specific parsing
abstract type AbstractMission end

# abstract supertype
abstract type AbstractMetadata end

missiontrait(T::Type{<:AbstractMetadata}) = error("No trait set for $(T)")
missiontrait(::M) where {M<:AbstractMetadata} = missiontrait(M)

# interface
trim_meta!(::AbstractMetadata, inds) = nothing

group_meta(::AbstractMetadata, mask, inds) = nothing

struct AncillaryResponse{T}
    bins_low::Vector{T}
    bins_high::Vector{T}
    spec_response::Vector{T}
end

# could be Response or Redistribution : how do we track this? 
struct ResponseMatrix{T}
    matrix::SparseMatrixCSC{T,Int}
    channels::Vector{Int}
    channel_bins_low::Vector{T}
    channel_bins_high::Vector{T}
    bins_low::Vector{T}
    bins_high::Vector{T}
end

# concrete type
mutable struct SpectralDataset{
    T,
    MetaType,
    PoissType,
    UnitType,
    AncType,
    BkgType,
    GroupType,
    VecType,
}
    # store high and low seperately
    # incase discontinuous dataset
    # - will there ever be discontinuous bins??
    bins_low::VecType
    bins_high::VecType

    _data::VecType
    _errors::VecType
    units::UnitType

    meta::MetaType
    poisson_errors::PoissType
    response::ResponseMatrix{T}
    ancillary::AncType
    background::BkgType

    channels::Vector{Int}
    grouping::GroupType
    quality::Vector{Int}

    mask::BitVector

    exposure_time::T
end

missiontrait(::SpectralDataset{T,M}) where {T,M} = missiontrait(M)

export SpectralDataset, ResponseMatrix
