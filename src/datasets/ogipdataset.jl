struct OGIPDataset{T,H} <: AbstractDataset
    data::SpectralData{T}
    paths::SpectralDataPaths
    observation_id::String
    exposure_id::String
    object::String
    header::H
end

function OGIPDataset(spec_path; T::Type = Float64, hdu = 2, kwargs...)
    paths = SpectralDataPaths(spec_path)
    config = StandardOGIPConfig(rmf_matrix_index = 2, rmf_energy_index = 3, T = T)

    header = read_fits_header(paths.spectrum; hdu = hdu)

    obs_id = haskey(header, "OBS_ID") ? header["OBS_ID"] : "[no observation id]"
    exposure_id = haskey(header, "EXP_ID") ? header["EXP_ID"] : "[no exposure id]"
    object = haskey(header, "OBJECT") ? header["OBJECT"] : "[no object]"

    data = SpectralData(paths, config; kwargs...)
    OGIPDataset(data, paths, obs_id, exposure_id, object, header)
end

@_forward_SpectralData_api OGIPDataset.data

function _printinfo(io, data::OGIPDataset{T}) where {T}
    descr = """OGIPDataset:
      . Object              : $(data.object)
      . Observation ID      : $(data.observation_id)
      . Exposure ID         : $(data.exposure_id)
    """
    print(io, descr)
    _printinfo(io, data.data)
end

export OGIPDataset
