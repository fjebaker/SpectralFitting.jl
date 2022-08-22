struct OGIP_Spectrum{T}
    exposure_time::T
    background_scale::T
    area_scale::T
    systematic_error::T
    poisson_error::Bool
    units::Symbol

    channels::Vector{Int}
    quality::Vector{Int}
    grouping::Vector{Int}
    # may be counts or rate depending
    values::Vector{T}
    stat_error::Vector{T}

    # metadata
    telescope::String
    instrument::String
end

function OGIP_Spectrum(fits::FITS, ::Type{T}) where {T}
    header = read_header(fits[2])
    is_poisson = header["POISSERR"] == "T" ? true : false
    instrument = header["INSTRUME"]
    telescope = header["TELESCOP"]
    exposure_time = T(header["EXPOSURE"])
    background_scale = T(header["BACKSCAL"])
    area_scale = T(header["AREASCAL"])
    sys_error = T(header["SYS_ERR"])

    channels = Int.(read(fits[2], "CHANNEL"))
    quality = Int.(read(fits[2], "QUALITY"))
    grouping = Int.(read(fits[2], "GROUPING"))

    column_names = FITSIO.colnames(fits[2])
    units, values = if "RATE" ∈ column_names
        :rate, T.(read(fits[2], "RATE"))
    else
        :counts, T.(read(fits[2], "COUNTS"))
    end
    stat_errors = if "STAT_ERR" ∈ column_names
        T.(read(fits[2], "STAT_ERR"))
    else
        T[]
    end

    OGIP_Spectrum(
        exposure_time,
        background_scale,
        area_scale,
        sys_error,
        is_poisson,
        units,
        channels,
        quality,
        grouping,
        values,
        stat_errors,
        telescope,
        instrument,
    )
end

struct OGIP_ARF{T}
    spec_response::Vector{T}
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
end

function OGIP_ARF(fits::FITS, ::Type{T}) where {T}
    OGIP_ARF(
        T.(read(fits[2], "SPECRESP")),
        T.(read(fits[2], "ENERG_LO")),
        T.(read(fits[2], "ENERG_HI")),
    )
end

struct OGIP_RMF_Channels{T}
    channels::Vector{Int}
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
end

function OGIP_RMF_Channels(fits::FITS, ::Type{T}) where {T}
    OGIP_RMF_Channels(
        Int.(read(fits[3], "CHANNEL")),
        T.(read(fits[3], "E_MIN")),
        T.(read(fits[3], "E_MAX")),
    )
end

struct OGIP_RMF_Matrix{T,M}
    Fchan::Matrix{Int}
    Nchan::Matrix{Int}
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
    matrix_rows::M
    first_channel::Int
    number_of_channels::Int
    number_of_energies::Int
end

function OGIP_RMF_Matrix(fits::FITS, first_channel, number_of_channels, ::Type{T}) where {T}
    e_lows = T.(read(fits[2], "ENERG_LO"))
    OGIP_RMF_Matrix(
        Int.(read(fits[2], "F_CHAN")),
        Int.(read(fits[2], "N_CHAN")),
        e_lows,
        T.(read(fits[2], "ENERG_HI")),
        read(fits[2], "MATRIX"),
        first_channel,
        number_of_channels,
        length(e_lows),
    )
end

function OGIP_RMF_Matrix(fits::FITS, T::Type)
    first_channel, number_of_channels = parse_rmf_fits_header(fits)
    OGIP_RMF_Matrix(fits, first_channel, number_of_channels, T)
end

function parse_rmf_fits_header(fits)
    header = read_header(fits[2])
    # what can we learn from the header?
    column_names = FITSIO.colnames(fits[2])
    f_chan_index = findfirst(==("F_CHAN"), column_names)

    if isnothing(f_chan_index)
        throw("F_CHAN column is missing from RMF.")
    end

    tl_min_key = "TLMIN$f_chan_index"
    channel_min = if haskey(header, tl_min_key)
        Int(header[tl_min_key])
    else
        @warn "No TLMIN key set in RMF header ($tl_min_key). Assuming channels start at 1."
        1
    end

    channel_num = if haskey(header, "DETCHANS")
        Int(header["DETCHANS"])
    else
        @warn "DETCHANS is not set in RMF header. Infering channel count from energy table."
        length(read(fits[3], "CHANNEL"))
    end

    channel_min, channel_num
end

struct OGIP_RMF{T}
    ogip_matrix::OGIP_RMF_Matrix{T}
    ogip_rmf_channels::OGIP_RMF_Channels{T}
end

OGIP_RMF(fits::FITS, T::Type) =
    OGIP_RMF(OGIP_RMF_Matrix(fits, T), OGIP_RMF_Channels(fits, T))

# utility constructors for path specs
function OGIP_Spectrum(path::String; T = Float64)
    fits = FITS(path)
    spec = OGIP_Spectrum(fits, T)
    close(fits)
    spec
end
function OGIP_ARF(path::String; T = Float64)
    fits = FITS(path)
    arf = OGIP_ARF(fits, T)
    close(fits)
    arf
end
function OGIP_RMF(path::String; T = Float64)
    fits = FITS(path)
    rmf = OGIP_RMF(fits, T)
    close(fits)
    rmf
end
