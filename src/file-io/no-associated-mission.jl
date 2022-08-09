
export drop_bad_channels!

mutable struct NoAssociatedMissionDatasetMeta{Q}
    quality::Q
    grp_path::String
    rm_path::String
end

missiontrait(::Type{<:NoAssociatedMissionDatasetMeta}) = NoAssociatedMission()

# create an alias
const SimpleSpectralDataset = SpectralDataset{<:NoAssociatedMissionDatasetMeta}

function trim_dataset_meta!(sdm::NoAssociatedMissionDatasetMeta, inds)
    sdm.quality = sdm.quality[inds]
end

drop_bad_channels!(sd::SimpleSpectralDataset) = trim_dataset!(sd, sd.meta.quality .== 0)

# must return these three things
function parse_grp_fits_file(::NoAssociatedMission, fits)
    counts = read(fits[2], "RATE")
    countserror = read(fits[2], "STAT_ERR")
    channels = read(fits[2], "CHANNEL")
    (counts, countserror, channels)
end

# must return these four things
function parse_rm_fits_file(::NoAssociatedMission, fits)
    rm_low_energy = read(fits[2], "ENERG_LO")
    rm_high_energy = read(fits[2], "ENERG_HI")
    F_chan = eachcol(read(fits[2], "F_CHAN"))
    N_chan = eachcol(read(fits[2], "N_CHAN"))
    matrix_lookup = read(fits[2], "MATRIX")

    low_energy_bins = read(fits[3], "E_MIN")
    high_energy_bins = read(fits[3], "E_MAX")
    channels = read(fits[3], "CHANNEL")

    matrix = build_matrix_response(
        rm_low_energy,
        rm_high_energy,
        F_chan,
        N_chan,
        matrix_lookup,
        low_energy_bins,
        high_energy_bins,
        channels,
    )

    (matrix, channels, low_energy_bins, high_energy_bins)
end

function meta_from_fits(::NoAssociatedMission, fits, rm_fits, grp_path, rm_path)
    qs = read(fits[2], "QUALITY")
    quality_vec::Vector{Int} = Int.(qs)
    NoAssociatedMissionDatasetMeta(quality_vec, grp_path, rm_path)
end
