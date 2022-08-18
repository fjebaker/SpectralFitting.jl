struct NoMission <: AbstractMission end

mutable struct BasicMetadata{Q} <: AbstractMetadata
    quality::Q
    grp_path::String
    rm_path::String
end

missiontrait(::Type{<:BasicMetadata}) = NoMission()

trim_meta!(sdm::BasicMetadata, inds) = sdm.quality = sdm.quality[inds]

function SpectralDataset(mission::NoMission, path, rm_path; T::Type = Float64)
    fits = FITS(path)
    fits_rm = FITS(rm_path)

    rm = ResponseMatrix(fits_rm, mission, T)

    qs = read(fits[2], "QUALITY")
    counts = read(fits[2], "RATE")
    countserror = read(fits[2], "STAT_ERR")
    channels = read(fits[2], "CHANNEL")

    quality_vec::Vector{Int} = Int.(qs)
    meta = BasicMetadata(quality_vec, path, rm_path)

    make_spectral_dataset(meta, rm, counts, countserror, channels, T)
end
