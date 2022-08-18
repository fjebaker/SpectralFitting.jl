
struct OGIP_RMF_Channels{T}
    channels::Vector{Int}
    energy_bins_low::Vector{T}
    energy_bins_high::Vector{T}
end

function OGIP_RMF_Channels(fits, ::Type{T}) where {T}
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

function OGIP_RMF_Matrix(fits, first_channel, number_of_channels, ::Type{T}) where {T}
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

function OGIP_RMF_Matrix(fits, T::Type)
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