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
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
    spec_response::Vector{T}
end

# could be Response or Redistribution : how do we track this? 
struct ResponseMatrix{T}
    matrix::SparseMatrixCSC{T,Int}
    channels::Vector{Int}
    channel_energy_bins_low::Vector{T}
    channel_energy_bins_high::Vector{T}
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
end

# concrete type
mutable struct SpectralDataset{T,M,P,K,A,B}
    # store high and low seperately
    # incase discontinuous dataset
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
    # counts and countserror have seperate type to support optional unit tracking
    counts::Vector{K}
    countserror::Vector{K}

    meta::M
    poisson_errors::P
    response::ResponseMatrix{T}
    ancillary::A
    background::B

    channels::Vector{Int}
    grouping::Vector{Int}
    quality::Vector{Int}

    mask::BitVector

    exposure_time::T
end

missiontrait(::SpectralDataset{T,M}) where {T,M} = missiontrait(M)

export SpectralDataset, ResponseMatrix
