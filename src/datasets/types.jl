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

# @enumx ErrorStatistics begin
#     Numeric
#     Poisson
#     Gaussian
#     Unknown
# end

struct SpectralFilePaths
    spectrum::Union{Missing,String}
    background::Union{Missing,String}
    response::Union{Missing,String}
    ancillary::Union{Missing,String}
end

path_assembler(::Type{NoMission}) = OGIP.read_paths_from_spectrum
path_assembler(::T) where {T<:AbstractMission} = path_assembler(T)

function SpectralFilePaths(; spectrum = "", background = "", response = "", ancillary = "")
    SpectralFilePaths(
        spectrum == "" ? missing : spectrum,
        background == "" ? missing : background,
        response == "" ? missing : response,
        ancillary == "" ? missing : ancillary,
    )
end

# struct AncillaryResponse{T}
#     bins_low::Vector{T}
#     bins_high::Vector{T}
#     effective_area::Vector{T}
# end

# TODO: could be Response or Redistribution : how do we track this? 
# struct ResponseMatrix{T}
#     matrix::SparseMatrixCSC{T,Int}
#     channels::Vector{Int}
#     channel_bins_low::Vector{T}
#     channel_bins_high::Vector{T}
#     bins_low::Vector{T}
#     bins_high::Vector{T}
# end

# struct Spectrum{T}
#     channels::Vector{Int}
#     quality::Vector{Int}
#     grouping::Vector{Int}

#     values::Vector{T}
#     unit_string::String

#     exposure_time::T
#     background_scale::T
#     area_scale::T

#     error_statistics::SpectralFitting.ErrorStatistics.T
#     errors::Union{Missing,Vector{T}}
#     systematic_error::T

#     telescope::String
#     instrument::String
# end


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

mutable struct SpectralDataset{UnitType,MetaType,T} <: AbstractDataset
    # store high and low seperately incase discontinuous dataset
    # TODO: will there ever be discontinuous bins??
    units::UnitType
    meta::MetaType

    bins_low::Vector{T}
    bins_high::Vector{T}

    spectrum::Spectrum{T}
    background::Union{Missing,Spectrum{T}}
    response::Union{Missing,ResponseMatrix{T}}
    ancillary::Union{Missing,AncillaryResponse{T}}

    mask::BitVector
end

# methods than can be subtypes for other dataset types
missiontrait(::SpectralDataset{U,M}) where {U,M} = missiontrait(M)

export SpectralDataset, ResponseMatrix, NoMission
