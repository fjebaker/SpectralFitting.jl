export AbstractXmmNewtonDevice, XmmNewton, XmmNewtonEPIC, XmmNewtonMeta

abstract type AbstractXmmNewtonDevice end
struct XmmEPIC <: AbstractXmmNewtonDevice end

struct XmmNewton{D} <: AbstractMission end
XmmNewton(::D) where {D<:AbstractXmmNewtonDevice} = XmmNewton{D}()

const XmmNewtonEPIC = XmmNewton{XmmEPIC}

struct XmmNewtonMeta{D} <: AbstractMetadata
    device::D
    paths::SpectralFilePaths
    observation_id::String
    exposure_id::String
    object::String
end

function observation_id(data::SpectralDataset{T,<:XmmNewtonMeta}) where {T}
    data.meta.observation_id
end
function observation_object(data::SpectralDataset{T,<:XmmNewtonMeta}) where {T}
    data.meta.object
end

SpectralFitting.missiontrait(::Type{<:XmmNewtonMeta{D}}) where {D} = XmmNewton(D())
SpectralFitting.path_assembler(::Type{<:XmmNewton}) = OGIP.read_paths_from_spectrum

function XmmNewtonMeta(device::AbstractXmmNewtonDevice, paths)
    fits = FITS(paths.spectrum)
    header = read_header(fits[1])
    close(fits)
    XmmNewtonMeta(
        device,
        paths,
        haskey(header, "OBS_ID") ? header["OBS_ID"] : "",
        haskey(header, "EXP_ID") ? header["EXP_ID"] : "",
        haskey(header, "OBJECT") ? header["OBJECT"] : "",
    )
end

function SpectralDataset(
    ::XmmNewton{D},
    paths::SpectralFilePaths,
    T::Type = Float64,
) where {D}
    config = StandardOGIPConfig(rmf_matrix_index = 2, rmf_energy_index = 3, T = T)
    meta = XmmNewtonMeta(D(), paths)

    dataset_from_ogip(paths, config, meta)
end
