# for future use: mission specific parsing
abstract type AbstractMission end
struct NoMission <: AbstractMission end

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
    effective_area::Vector{T}
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

abstract type AbstractDataset end
missiontrait(::AbstractDataset) = NoMission()
target_vector(data::AbstractDataset) = error(
    "No implemented for $(typeof(data)) yet. This function should return the y vector for the fit.",
)
target_variance(data::AbstractDataset) = error(
    "No implemented for $(typeof(data)) yet. This function should return the variance on the y vector for the fit.",
)
domain_vector(data::AbstractDataset) = error(
    "No implemented for $(typeof(data)) yet. This function should return the x vector for the fit.",
)
function _lazy_folded_invokemodel(model::AbstractSpectralModel, ::AbstractDataset)
    # no response data here
    (x, params) -> invokemodel(x, model, params)
end

mutable struct SpectralDataset{
    T,
    MetaType,
    PoissType,
    UnitType,
    AncType,
    BkgType,
    GroupType,
    VecType,
} <: AbstractDataset
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

# methods than can be subtypes for other dataset types
missiontrait(::SpectralDataset{T,M}) where {T,M} = missiontrait(M)
target_vector(data::SpectralDataset) = data.rate
target_variance(data::SpectralDataset) = data.rateerror .^ 2
domain_vector(data::SpectralDataset) = domain_vector(data.response)
function _lazy_folded_invokemodel(model::AbstractSpectralModel, data::SpectralDataset)
    ΔE = data.energy_bin_widths
    # pre-mask the response matrix to ensure channel out corresponds to the active data points
    R = fold_ancillary(data)[data.mask, :]
    # pre-allocate the output 
    wrapped = (energy, params) -> begin
        flux = invokemodel(energy, model, params)
        flux = (R * flux)
        @. flux = flux / ΔE
    end
    wrapped
end

export SpectralDataset, ResponseMatrix, NoMission
