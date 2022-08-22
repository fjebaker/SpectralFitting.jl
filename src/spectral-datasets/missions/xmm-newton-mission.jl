export AbstractXmmNewtonDevice, EPIC, XmmNewton, XmmNewtonEPIC, XmmNewtonMeta

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
end

SpectralFitting.missiontrait(::Type{<:XmmNewtonMeta{D}}) where {D} = XmmNewton(D())

function SpectralDataset(
    ::XmmNewton{D},
    path,
    rm_path,
    arf_path;
    T::Type = Float64,
) where {D}
    spec = OGIP_Spectrum(path; T)
    rm = OGIP_RMF(rm_path; T)
    arf = OGIP_ARF(arf_path; T)
    meta = XmmNewtonMeta(D(), path, rm_path, arf_path)
    SpectralDataset(spec, rm, arf, nothing, meta)
end
