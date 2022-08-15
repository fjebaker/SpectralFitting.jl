mutable struct NoAssociatedMissionMeta{Q} <: AbstractSpectralDatasetMeta
    quality::Q
    grp_path::String
    rm_path::String
end

missiontrait(::Type{<:NoAssociatedMissionMeta}) = NoAssociatedMission()

trim_dataset_meta!(sdm::NoAssociatedMissionMeta, inds) =
    sdm.quality = sdm.quality[inds]

drop_bad_channels!(sd, ::NoAssociatedMission) = trim_dataset!(sd, sd.meta.quality .== 0)

function load_spectral_dataset(mission::NoAssociatedMission, path, rm_path; T::Type = Float64)
    fits = FITS(path)
    fits_rm = FITS(rm_path)

    rm = load_response_matrix(fits_rm, mission, T)

    qs = read(fits[2], "QUALITY")
    counts = read(fits[2], "RATE")
    countserror = read(fits[2], "STAT_ERR")
    channels = read(fits[2], "CHANNEL")

    quality_vec::Vector{Int} = Int.(qs)
    meta = NoAssociatedMissionMeta(quality_vec, path, rm_path)

    make_spectral_dataset(meta, rm, counts, countserror, channels, T)
end