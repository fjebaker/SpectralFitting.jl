module OGIP
import SpectralFitting
using FITSIO

export AbstractOGIPConfig, StandardOGIPConfig

abstract type AbstractOGIPConfig{T} end
rmf_matrix_index(c::AbstractOGIPConfig) = error("Not implemented for $(typeof(c))")
rmf_energy_index(c::AbstractOGIPConfig) = error("Not implemented for $(typeof(c))")

struct StandardOGIPConfig{T} <: AbstractOGIPConfig{T}
    rmf_matrix_index::Int
    rmf_energy_index::Int
    function StandardOGIPConfig(;
        rmf_matrix_index = 3,
        rmf_energy_index = 2,
        T::Type = Float64,
    )
        @assert rmf_matrix_index != rmf_energy_index
        new{T}(rmf_matrix_index, rmf_energy_index)
    end
end
rmf_matrix_index(c::StandardOGIPConfig) = c.rmf_matrix_index
rmf_energy_index(c::StandardOGIPConfig) = c.rmf_energy_index

struct MissingHeader <: Exception
    header::String
end
Base.showerror(io::IO, e::MissingHeader) = print(io, "Header: '$(e.header)' is not defined")

# structs
struct RMFHeader
    first_channel::Int
    num_channels::Int
end

struct RMFMatrix{V,T,M}
    f_chan::V
    n_chan::V
    bins_low::Vector{T}
    bins_high::Vector{T}
    matrix::M
    header::RMFHeader
end

struct RMFChannels{T}
    channels::Vector{Int}
    bins_low::Vector{T}
    bins_high::Vector{T}
end

# functions
function parse_rmf_header(table::TableHDU)
    header = FITSIO.read_header(table)
    columns = FITSIO.colnames(table)

    findex = findfirst(==("F_CHAN"), columns)
    if isnothing(findex)
        throw(MissingHeader("F_CHAIN"))
    end

    tlindex = "TLMIN$findex"
    first_channel = if haskey(header, tlindex)
        Int(header[tlindex])
    else
        @warn "No TLMIN key set in RMF header ($tl_min_key). Assuming channels start at 1."
        1
    end
    num_channels = if haskey(header, "DETCHANS")
        Int(header["DETCHANS"])
    else
        @warn "DETCHANS is not set in RMF header. Infering channel count from table length."
        -1
    end
    RMFHeader(first_channel, num_channels)
end

function read_rmf_channels(table::TableHDU, T::Type)
    channels = Int.(read(table, "CHANNEL"))
    energy_low = T.(read(table, "E_MIN"))
    energy_high = T.(read(table, "E_MAX"))
    RMFChannels(channels, energy_low, energy_high)
end

function read_rmf_matrix(table::TableHDU, header::RMFHeader, T::Type)
    energy_low = convert.(T, read(table, "ENERG_LO"))
    energy_high = convert.(T, read(table, "ENERG_HI"))
    f_chan_raw = read(table, "F_CHAN")
    n_chan_raw = read(table, "N_CHAN")
    matrix_raw = read(table, "MATRIX")
    RMFMatrix(f_chan_raw, n_chan_raw, energy_low, energy_high, matrix_raw, header)
end

function read_rmf(path::String, ogip_config::AbstractOGIPConfig; T::Type = Float64)
    fits = FITS(path)

    i = rmf_matrix_index(ogip_config)
    j = rmf_energy_index(ogip_config)

    (header, res, channels) = try
        header = parse_rmf_header(fits[i])
        res = read_rmf_matrix(fits[i], header, T)
        channels::RMFChannels{T} = read_rmf_channels(fits[j], T)
        (header, res, channels)
    catch e
        throw(e)
        # ... 
    finally
        close(fits)
    end

    channels
end

function SpectralFitting.ResponseMatrix(rmf::RMFMatrix, channels::RMFChannels) end

end # module

using .OGIP
export StandardOGIPConfig
