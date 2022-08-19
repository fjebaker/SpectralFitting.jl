export AbstractXmmNewtonDevice, EPIC, XmmNewton, XmmNewtonEPIC, XmmNewtonMeta

struct AncillaryResponse{T}
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
    spec_response::Vector{T}
end

abstract type AbstractXmmNewtonDevice end
struct XmmEPIC <: AbstractXmmNewtonDevice end

struct XmmNewton{D} <: AbstractMission end
XmmNewton(::D) where {D<:AbstractXmmNewtonDevice} = XmmNewton{D}()

const XmmNewtonEPIC = XmmNewton{XmmEPIC}

mutable struct XmmNewtonMeta{D,A} <: AbstractMetadata
    device::D
    quality::Vector{Int}
    grouping::Vector{Int}
    ancillary_response::A
    grp_path::String
    rm_path::String
    arf_path::String
    exposure_time::Float64
end

SpectralFitting.missiontrait(::Type{<:XmmNewtonMeta{D}}) where {D} = XmmNewton(D())

function trim_meta!(sdm::XmmNewtonMeta, inds)
    sdm.quality = sdm.quality[inds]
    sdm.grouping = sdm.grouping[inds]
end

function group_meta(sdm::XmmNewtonMeta, mask, inds)
    N = length(inds) - 1
    new_quality = group_quality_vector(mask, sdm.quality, inds)
    XmmNewtonMeta(
        sdm.device,
        new_quality,
        # update grouping
        ones(Int, N),
        sdm.ancillary_response,
        sdm.grp_path,
        sdm.rm_path,
        sdm.arf_path,
        sdm.exposure_time,
    )
end

function SpectralDataset(
    mission::XmmNewton{D},
    path,
    rm_path,
    arf_path;
    T::Type = Float64
) where {D}
    fits = FITS(path)
    fits_rm = FITS(rm_path)

    rm = ResponseMatrix(fits_rm, mission, T)

    qs = read(fits[2], "QUALITY")

    exp_time = read_header(fits[2])["EXPOSURE"]

    counts = read(fits[2], "COUNTS")

    countserror = sqrt.(counts)
    channels = read(fits[2], "CHANNEL")
    grouping = read(fits[2], "GROUPING")

    close(fits)
    close(fits_rm)

    quality_vec::Vector{Int} = Int.(qs)
    grouping_vec::Vector{Int} = Int.(grouping)
    meta =
        XmmNewtonMeta(D(), quality_vec, grouping_vec, (), path, rm_path, arf_path, exp_time)

    SpectralDataset(meta, rm, counts, countserror, channels; T = T)
end
