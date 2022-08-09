export load_spectral_dataset

function load_spectral_dataset(
    grp_path,
    rm_path;
    mission::AbstractMissionTrait = NoAssociatedMission(),
    T::Type = Float64,
)
    f_grp = FITS(grp_path)
    f_rm = FITS(rm_path)

    (counts, countserror, channels) = parse_grp_fits_file(mission, f_grp)
    (matrix, rm_channels, rm_low_energy_bins, rm_high_energy_bins) =
        parse_rm_fits_file(mission, f_rm)

    meta = meta_from_fits(mission, f_grp, f_rm, grp_path, rm_path)

    close(f_grp)
    close(f_rm)

    rm::ResponseMatrix{T} = ResponseMatrix(
        matrix,
        Int.(rm_channels),
        T.(rm_low_energy_bins),
        T.(rm_high_energy_bins),
    )

    low_energy_bins, high_energy_bins = augmented_energy_channels(channels, rm)
    counts_vec::Vector{T} = T.(counts)
    countserr_vec::Vector{T} = T.(countserror)
    channels_vec::Vector{Int} = Int.(channels)

    SpectralDataset(
        meta,
        rm,
        low_energy_bins,
        high_energy_bins,
        counts_vec,
        countserr_vec,
        channels_vec,
        false,
    )
end
