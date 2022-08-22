struct NoMission <: AbstractMission end

struct BasicMetadata{Q} <: AbstractMetadata
    grp_path::String
    rm_path::String
    telescope::String
    instrument::String
end

missiontrait(::Type{<:BasicMetadata}) = NoMission()

function SpectralDataset(::NoMission, path, rm_path; T::Type = Float64)
    spec = OGIP_Spectrum(path; T)
    rmf = OGIP_RMF(rm_path; T)
    meta = BasicMetadata(path, rm_path, spec.telescope, spec.instrument)
    SpectralDataset(spec, rmf, nothing, nothing, meta)
end
