struct OGIPDataset{T,H} <: AbstractDataset
    data::SpectralData{T}
    paths::SpectralDataPaths
    observation_id::String
    exposure_id::String
    object::String
    header::H
end

function load_ogip_dataset(
    spec_path;
    hdu = 2,
    background = missing,
    response = missing,
    ancillary = missing,
    kwargs...,
)
    paths = SpectralDataPaths(
        spec_path;
        background = background,
        response = response,
        ancillary = ancillary,
    )
    config = StandardOGIPConfig(; kwargs...)

    header = read_fits_header(paths.spectrum; hdu = hdu)

    obs_id = haskey(header, "OBS_ID") ? header["OBS_ID"] : "[no observation id]"
    exposure_id = haskey(header, "EXP_ID") ? header["EXP_ID"] : "[no exposure id]"
    object = haskey(header, "OBJECT") ? header["OBJECT"] : "[no object]"

    data = SpectralData(paths, config)
    (data, paths, obs_id, exposure_id, object, header)
end

OGIPDataset(spec_path; kwargs...) = OGIPDataset(load_ogip_dataset(spec_path; kwargs...)...)

@_forward_SpectralData_api OGIPDataset.data

make_label(d::OGIPDataset) = d.observation_id

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
