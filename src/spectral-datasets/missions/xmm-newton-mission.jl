export AbstractXmmNewtonDevice, XmmNewton, XmmNewtonEPIC, XmmNewtonMeta

abstract type AbstractXmmNewtonDevice end
struct XmmEPIC <: AbstractXmmNewtonDevice end

struct XmmNewton{D} <: AbstractMission end
XmmNewton(::D) where {D<:AbstractXmmNewtonDevice} = XmmNewton{D}()

const XmmNewtonEPIC = XmmNewton{XmmEPIC}

struct XmmNewtonMeta{D} <: AbstractMetadata
    device::D
    grp_path::String
    rm_path::String
    arf_path::String
    observation_id::String
    exposure_id::String
    object::String
end

function observation_id(data::SpectralDataset{T,<:XmmNewtonMeta}) where {T}
    data.meta.observation_id
end

SpectralFitting.missiontrait(::Type{<:XmmNewtonMeta{D}}) where {D} = XmmNewton(D())

function XmmNewtonMeta(device::AbstractXmmNewtonDevice, path, rm_path, arf_path)
    fits = FITS(path)
    header = read_header(fits[1])
    close(fits)
    XmmNewtonMeta(
        device,
        path,
        rm_path,
        arf_path,
        haskey(header, "OBS_ID") ? header["OBS_ID"] : "",
        haskey(header, "EXP_ID") ? header["EXP_ID"] : "",
        haskey(header, "OBJECT") ? header["OBJECT"] : "",
    )
end

function SpectralDataset(
    ::XmmNewton{D},
    path,
    rm_path,
    arf_path
    ;
    T::Type = Float64,
) where {D}
    spec = OGIP_GroupedEvents(path; T)
    rm = OGIP_RMF(rm_path; T)
    arf = OGIP_ARF(arf_path; T)
    meta = XmmNewtonMeta(D(), path, rm_path, arf_path)
    SpectralDataset(SpectralUnits.u"counts", spec, rm, arf, nothing, meta)
end
