
struct NuStarData{T,H} <: AbstractDataset
    data::SpectralData{T}
    paths::SpectralDataPaths
    observation_id::String
    exposure_id::String
    object::String
    header::H
end

NuStarData(spec_path; rmf_matrix_index = 3, rmf_energy_index = 2, kwargs...) = NuStarData(
    load_ogip_dataset(
        spec_path;
        rmf_matrix_index = rmf_matrix_index,
        rmf_energy_index = rmf_energy_index,
        kwargs...,
    )...,
)

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


abstract type AbstractXmmNewtonDevice end
struct XmmEPIC <: AbstractXmmNewtonDevice end

struct XmmData{T,H,D} <: AbstractDataset
    device::D
    data::SpectralData{T}
    paths::SpectralDataPaths
    observation_id::String
    exposure_id::String
    object::String
    header::H
end

XmmData(
    device::AbstractXmmNewtonDevice,
    spec_path;
    rmf_matrix_index = 2,
    rmf_energy_index = 3,
    kwargs...,
) = XmmData(
    device,
    load_ogip_dataset(
        spec_path;
        rmf_matrix_index = rmf_matrix_index,
        rmf_energy_index = rmf_energy_index,
        kwargs...,
    )...,
)

make_label(data::XmmData) = data.observation_id

@_forward_SpectralData_api XmmData.data

function Base.show(io::IO, @nospecialize(data::XmmData{T})) where {T}
    print(io, "XmmData[dev=$(data.device),obs_id=$(data.observation_id)]")
end

function _printinfo(io, data::XmmData{T}) where {T}
    descr = """XmmData for $(Base.typename(typeof(data.device)).name):
      . Object              : $(data.object)
      . Observation ID      : $(data.observation_id)
      . Exposure ID         : $(data.exposure_id)
    """
    print(io, descr)
    _printinfo(io, data.data)
end


export NuStarData, XmmData, XmmEPIC
