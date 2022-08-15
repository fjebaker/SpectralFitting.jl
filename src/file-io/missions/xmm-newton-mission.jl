export AbstractXMMNewtonDevice, EPIC, XMMNewton, XMMNewtonEPIC, XMMNewtonMeta

struct AncillaryResponse{T}
    low_energy_bins::Vector{T}
    high_energy_bins::Vector{T}
    spec_response::Vector{T}
end

abstract type AbstractXMMNewtonDevice end
struct EPIC <: AbstractXMMNewtonDevice end

struct XMMNewton{D} <: AbstractMissionTrait end
XMMNewton(::D) where {D<:AbstractXMMNewtonDevice} = XMMNewton{D}()

const XMMNewtonEPIC = XMMNewton{EPIC}

mutable struct XMMNewtonMeta{D,Q,G,A} <: AbstractSpectralDatasetMeta
    device::D
    quality::Q
    grouping::G
    ancillary_response::A
    grp_path::String
    rm_path::String
    arf_path::String
    exposure_time::Float64
    grouped::Bool
end

SpectralFitting.missiontrait(::Type{<:XMMNewtonMeta{D}}) where {D} = D()

is_grouped(sd::SpectralDataset{<:XMMNewtonMeta}) = sd.meta.grouped

function trim_dataset_meta!(sdm::XMMNewtonMeta, inds)
    sdm.quality = sdm.quality[inds]
    sdm.grouping = sdm.grouping[inds]
end

function group_dataset_meta!(sdm::XMMNewtonMeta, inds, T)
    sdm.grouped = true
    
    N = length(inds) - 1

    new_quality = zeros(T, N)
    grouping_indices_callback(inds) do (i, index1, index2)
        new_quality[i] = @views round(Int, Statistics.mean(sdm.quality[index1:index2]))
    end

    sdm.quality = new_quality
end

drop_bad_channels!(sd, ::XMMNewton) = trim_dataset!(sd, sd.meta.quality .== 0)

function load_spectral_dataset(mission::XMMNewton, path, rm_path, arf_path; T::Type = Float64)
    fits = FITS(path)
    fits_rm = FITS(rm_path)

    rm = load_response_matrix(fits_rm, mission, T)

    qs = read(fits[2], "QUALITY")

    exp_time = read_header(fits[2])["EXPOSURE"]

    # turn into counts per second
    counts = read(fits[2], "COUNTS") ./ exp_time
    
    countserror = sqrt.(counts)
    channels = read(fits[2], "CHANNEL")
    grouping = read(fits[2], "GROUPING")

    quality_vec::Vector{Int} = Int.(qs)
    meta = XMMNewtonMeta(mission, quality_vec, grouping, (), path, rm_path, arf_path, exp_time, false)

    make_spectral_dataset(meta, rm, counts, countserror, channels, T)
end