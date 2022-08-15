export load_spectral_dataset

function load_response_matrix(fits, mission, T)
    (matrix, rm_channels, rm_low_energy_bins, rm_high_energy_bins) =
        parse_rm_fits_file(mission, fits, T)

    rm::ResponseMatrix{T} = ResponseMatrix(
        matrix,
        Int.(rm_channels),
        T.(rm_low_energy_bins),
        T.(rm_high_energy_bins),
    )

    return rm
end

function make_spectral_dataset(meta::AbstractSpectralDatasetMeta, rm::ResponseMatrix, counts, countserror, channels, T)
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

load_spectral_dataset(::AbstractMissionTrait, args...; kwargs...) = error("First argument must be a mission trait.")
