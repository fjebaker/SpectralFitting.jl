export NoMission

struct NoMission <: AbstractMission end

struct BasicMetadata <: AbstractMetadata
    grp_path::String
    rm_path::String
    telescope::String
    instrument::String
end

missiontrait(::Type{<:BasicMetadata}) = NoMission()
function observation_id(data::SpectralDataset{T,<:BasicMetadata}) where {T}
    "noID"
end
function observation_object(data::SpectralDataset{T,<:BasicMetadata}) where {T}
    "Unkown"
end

function SpectralDataset(::NoMission, path, rm_path; T::Type = Float64, units = nothing)
    spec = OGIP_GroupedEvents(path; T)
    rmf = OGIP_RMF(rm_path; T)
    meta = BasicMetadata(path, rm_path, spec.telescope, spec.instrument)
    _units = if isnothing(units) # if no units given, have to guess them from the data
        infer_units(spec.units)
    else
        units
    end
    SpectralDataset(_units, spec, rmf, nothing, nothing, meta)
end
