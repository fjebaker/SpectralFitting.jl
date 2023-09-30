
struct NuStarData{T} <: AbstractDataset
    data::SpectralData{T}
    paths::SpectralDataPaths
    observation_id::String
    exposure_id::String
    object::String
end

function NuStarData(spec_path; T::Type = Float64, kwargs...)
    paths = SpectralDataPaths(spec_path)
    config = StandardOGIPConfig(rmf_matrix_index = 3, rmf_energy_index = 2, T = T)

    # read metadata
    fits = FITS(paths.spectrum)
    header = read_header(fits[1])
    close(fits)

    obs_id = haskey(header, "OBS_ID") ? header["OBS_ID"] : "[no observation id]"
    exposure_id = haskey(header, "EXP_ID") ? header["EXP_ID"] : "[no exposure id]"
    object = haskey(header, "OBJECT") ? header["OBJECT"] : "[no object]"

    data = SpectralData(paths, config; kwargs...)
    NuStarData(data, paths, obs_id, exposure_id, object)
end

make_label(data::NuStarData) = data.observation_id

@_forward_SpectralData_api NuStarData.data

function Base.show(io::IO, @nospecialize(data::NuStarData{T})) where {T}
    print(io, "NuStarData[obs_id=$(data.observation_id)]")
end

function _printinfo(io, data::NuStarData{T}) where {T}
    descr = """NuStarData:
      . Object              : $(data.object)
      . Observation ID      : $(data.observation_id)
      . Exposure ID         : $(data.exposure_id)
    """
    print(io, descr)
    _printinfo(io, data.data)
end


export NuStarData
