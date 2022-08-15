export AbstractMissionTrait,
    NoAssociatedMission,
    missiontrait,
    drop_bad_channels!

# for future use: mission specific parsing
abstract type AbstractMissionTrait end
struct NoAssociatedMission <: AbstractMissionTrait end

missiontrait(T::Type{<:AbstractSpectralDatasetMeta}) = error("No trait set for $(T)")
missiontrait(::M) where {M<:AbstractSpectralDatasetMeta} = missiontrait(M)
missiontrait(::SpectralDataset{M}) where {M} = missiontrait(M)

# interface depending on the meta
drop_bad_channels!(sd, M::AbstractMissionTrait; kwargs...) =
    nothing
drop_bad_channels!(sd::SpectralDataset; kwargs...) = drop_bad_channels!(sd, missiontrait(sd); kwargs...)

# must return these four things
function parse_rm_fits_file(mission::AbstractMissionTrait, fits, T)
    # from matrix table
    rm_low_energy = read(fits[2], "ENERG_LO")
    rm_high_energy = read(fits[2], "ENERG_HI")
    F_chan = eachcol(read(fits[2], "F_CHAN"))
    N_chan = eachcol(read(fits[2], "N_CHAN"))
    matrix_lookup = read(fits[2], "MATRIX")
    # from energy table
    low_energy_bins = read(fits[3], "E_MIN")
    high_energy_bins = read(fits[3], "E_MAX")
    channels = read(fits[3], "CHANNEL")

    N = length(channels)

    offset = mission == NoAssociatedMission() ? 1 : 0
    matrix = spzeros(T, N, N)

    build_matrix_response!(
        matrix,
        # from rm table
        rm_low_energy,
        rm_high_energy,
        F_chan,
        N_chan,
        matrix_lookup,
        # from energy table
        low_energy_bins,
        high_energy_bins,
        channels,
        offset
    )

    (matrix, channels, low_energy_bins, high_energy_bins)
end

# include missions
include("no-associated-mission.jl")
include("xmm-newton-mission.jl")