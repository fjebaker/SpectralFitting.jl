export NoMission

struct NoMission <: AbstractMission end

struct BasicMetadata <: AbstractMetadata
    grp_path::String
    rm_path::String
    telescope::String
    instrument::String
end

missiontrait(::Type{<:BasicMetadata}) = NoMission()

function SpectralDataset(::NoMission, path, rm_path; T::Type = Float64, units = nothing)
    spec = OGIP_Events(path; T)
    rmf = OGIP_RMF(rm_path; T)
    meta = BasicMetadata(path, rm_path, spec.telescope, spec.instrument)
    _units = if isnothing(units) # if no units given, have to guess them from the data
        infer_units(spec.units)
    else
        units
    end
    SpectralDataset(_units, spec, rmf, nothing, nothing, meta)
end
