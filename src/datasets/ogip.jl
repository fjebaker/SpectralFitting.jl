module OGIP

import SpectralFitting
import SpectralFitting: SpectralUnits
using FITSIO
using SparseArrays

export AbstractOGIPConfig, StandardOGIPConfig

abstract type AbstractOGIPConfig{T} end
rmf_matrix_index(c::AbstractOGIPConfig) = error("Not implemented for $(typeof(c))")
rmf_energy_index(c::AbstractOGIPConfig) = error("Not implemented for $(typeof(c))")

struct StandardOGIPConfig{T} <: AbstractOGIPConfig{T}
    rmf_matrix_index::Int
    rmf_energy_index::Int
    rmf_matrix_key::String
    function StandardOGIPConfig(;
        rmf_matrix_index = 3,
        rmf_energy_index = 2,
        rmf_matrix_key = "SPECRESP",
        T::Type = Float64,
    )
        @assert rmf_matrix_index != rmf_energy_index
        new{T}(rmf_matrix_index, rmf_energy_index, rmf_matrix_key)
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

function _parse_any(::Type{T}, @nospecialize(value::V))::T where {T,V}
    if V <: AbstractString
        parse(T, value)
    else
        convert(T, value)
    end
end

function _string_boolean(@nospecialize(value::V))::Bool where {V}
    if V <: AbstractString
        if value == "F"
            false
        elseif value == "T"
            true
        else
            @warn("Unknown boolean string: $(value)")
            false
        end
    else
        value
    end
end

# functions
function parse_rmf_header(table::TableHDU)
    header = FITSIO.read_header(table)
    columns = FITSIO.colnames(table)

    findex = findfirst(==("F_CHAN"), columns)
    if isnothing(findex)
        throw(MissingHeader("F_CHAN"))
    end

    tlindex = "TLMIN$findex"
    first_channel = if haskey(header, tlindex)
        _parse_any(Int, header[tlindex])
    else
        @warn "No TLMIN key set in RMF header ($tl_min_key). Assuming channels start at 1."
        1
    end
    num_channels = if haskey(header, "DETCHANS")
        _parse_any(Int, header["DETCHANS"])
    else
        @warn "DETCHANS is not set in RMF header. Infering channel count from table length."
        -1
    end
    RMFHeader(first_channel, num_channels)
end

function read_rmf_channels(table::TableHDU, T::Type)
    channels = _parse_any.(Int, read(table, "CHANNEL"))
    energy_low = _parse_any.(T, read(table, "E_MIN"))
    energy_high = _parse_any.(T, read(table, "E_MAX"))
    RMFChannels(channels, energy_low, energy_high)
end

function _chan_to_vectors(chan::Matrix)
    map(eachcol(chan)) do column
        i = findfirst(==(0), column)
        # if no zeroes, return full column
        if isnothing(i)
            column
        else
            # exclude the zero from our selection
            i = i > 1 ? i - 1 : i
            column[1:i]
        end
    end
end

function _translate_channel_array(channel)
    if channel isa Matrix
        _chan_to_vectors(channel)
    elseif eltype(channel) <: AbstractVector
        channel
    else
        map(i -> [i], channel)
    end
end

function _adapt_matrix_type(T::Type, mat::M) where {M}
    if eltype(M) <: AbstractVector
        map(row -> convert.(T, row), mat)
    elseif M <: AbstractMatrix
        map(row -> convert.(T, row), eachcol(mat))
    end
end

function read_rmf_matrix(table::TableHDU, header::RMFHeader, T::Type)
    energy_low = convert.(T, read(table, "ENERG_LO"))
    energy_high = convert.(T, read(table, "ENERG_HI"))
    f_chan_raw = read(table, "F_CHAN")
    n_chan_raw = read(table, "N_CHAN")
    matrix_raw = read(table, "MATRIX")

    # type stable: convert to common vector of vector format
    f_chan::Vector{Vector{Int}} = _translate_channel_array(f_chan_raw)
    n_chan::Vector{Vector{Int}} = _translate_channel_array(n_chan_raw)

    @show typeof(matrix_raw)
    @show size(f_chan), size(matrix_raw)

    RMFMatrix(
        f_chan,
        n_chan,
        energy_low,
        energy_high,
        _adapt_matrix_type(T, matrix_raw),
        header,
    )
end

function _read_fits_and_close(f, path)
    fits = FITS(path)
    res = try
        f(fits)
    catch e
        throw(e)
    finally
        close(fits)
    end
    res
end

function read_rmf(path::String, ogip_config::AbstractOGIPConfig{T}) where {T}
    i = rmf_matrix_index(ogip_config)
    j = rmf_energy_index(ogip_config)

    (header, rmf, channels::RMFChannels{T}) = _read_fits_and_close(path) do fits
        hdr = parse_rmf_header(fits[i])
        _rmf = read_rmf_matrix(fits[i], hdr, T)
        _channels = read_rmf_channels(fits[j], T)
        (hdr, _rmf, _channels)
    end

    _build_reponse_matrix(header, rmf, channels, T)
end

function read_ancillary_response(path::String, ogip_config::AbstractOGIPConfig{T}) where {T}
    fits = FITS(path)
    (bins_low, bins_high, effective_area) = _read_fits_and_close(path) do fits
        area::Vector{T} = convert.(T, read(fits[2], ogip_config.rmf_matrix_key))
        lo::Vector{T} = convert.(T, read(fits[2], "ENERG_LO"))
        hi::Vector{T} = convert.(T, read(fits[2], "ENERG_HI"))
        (lo, hi, area)
    end
    SpectralFitting.AncillaryResponse(bins_low, bins_high, effective_area)
end

function _build_reponse_matrix(
    header::RMFHeader,
    rmf::RMFMatrix,
    channels::RMFChannels,
    T::Type,
)
    R = spzeros(T, header.num_channels, length(rmf.bins_low))
    build_response_matrix!(R, rmf.f_chan, rmf.n_chan, rmf.matrix, header.first_channel)
    SpectralFitting.ResponseMatrix(
        R,
        channels.channels,
        channels.bins_low,
        channels.bins_high,
        rmf.bins_low,
        rmf.bins_high,
    )
end

function build_response_matrix!(
    R,
    f_chan::Vector,
    n_chan::Vector,
    matrix_rows::Vector,
    first_channel,
)
    for (i, (F, N)) in enumerate(zip(f_chan, n_chan))
        M = matrix_rows[i]
        index = 1
        for (first, len) in zip(F, N)
            if len == 0
                break
            end
            first -= first_channel
            @views R[first+1:first+len, i] .= M[index:index+len-1]
            index += len
        end
    end
end

function _read_exposure_time(header)
    if get(header, "EXPOSURE", 0.0) != 0.0
        return header["EXPOSURE"]
    end
    if get(header, "TELAPSE", 0.0) != 0.0
        return header["TELAPSE"]
    end
    # maybe time stops given
    if (get(header, "TSTART", 0.0) != 0.0) && (get(header, "TSTOP", 0.0) != 0.0)
        return header["TSTOP"] - header["TSTART"]
    end
    @warn "Cannot find or infer exposure time."
    0.0
end

function _get_stable(::Type{T}, header, name, default)::T where {T}
    get(header, name, T(default))
end

function read_spectrum(path, config::AbstractOGIPConfig{T}) where {T}
    info = _read_fits_and_close(path) do fits
        header = read_header(fits[2])
        # if not set, assume not poisson errors
        is_poisson = _string_boolean(get(header, "POISSERR", false))
        # read general infos
        instrument = header["INSTRUME"]
        telescope = header["TELESCOP"]
        exposure_time = T(_read_exposure_time(header))
        background_scale = _get_stable(T, header, "BACKSCAL", one(T))
        area_scale = _get_stable(T, header, "AREASCAL", one(T))
        sys_error = _get_stable(T, header, "SYS_ERR", zero(T))

        column_names = FITSIO.colnames(fits[2])

        channels::Vector{Int} = convert.(Int, read(fits[2], "CHANNEL"))

        quality::Vector{Int} = if "QUALITY" ∈ column_names
            convert.(Int, read(fits[2], "QUALITY"))
        else
            zeros(Int, size(channels))
        end

        grouping::Vector{Int} = if "GROUPING" ∈ column_names
            convert.(Int, read(fits[2], "GROUPING"))
        else
            ones(Int, size(channels))
        end

        units::SpectralUnits.RateOrCount, values::Vector{T} = if "RATE" ∈ column_names
            SpectralUnits._rate(), convert.(T, read(fits[2], "RATE"))
        else
            SpectralUnits._counts(), convert.(T, read(fits[2], "COUNTS"))
        end

        stat, _errors = if "STAT_ERR" ∈ column_names
            if is_poisson
                @warn "Both STAT_ERR column present and POISSERR flag set. Using STAT_ERR."
            end
            SpectralFitting.ErrorStatistics.Numeric, convert.(T, read(fits[2], "STAT_ERR"))
        elseif is_poisson
            SpectralFitting.ErrorStatistics.Poisson,
            @. T(SpectralFitting.count_error(values, 1.0))
        else
            @warn "Unknown error statistics. Setting zero for all."
            SpectralFitting.ErrorStatistics.Unknown, T[0 for _ in values]
        end

        SpectralFitting.Spectrum(
            channels,
            quality,
            grouping,
            values,
            units,
            exposure_time,
            background_scale,
            area_scale,
            stat,
            _errors,
            sys_error,
            telescope,
            instrument,
        )
    end
    info
end

function read_background(path::String, config::AbstractOGIPConfig{T}) where {T}
    read_spectrum(path, config)
end

function read_paths_from_spectrum(path::String)
    header = _read_fits_and_close(path) do fits
        read_header(fits[2])
    end
    # extract path information from header

    possible_ext = splitext(path)[2]
    response_path = read_filename(header, "RESPFILE", path, ".rmf", ".rsp")
    ancillary_path = read_filename(header, "ANCRFILE", path, possible_ext)
    background_path = read_filename(header, "BACKFILE", path, possible_ext)
    (background_path, response_path, ancillary_path)
end

function read_filename(header, entry, parent, exts...)
    data_directory = Base.dirname(parent)
    parent_name = basename(parent)
    if haskey(header, entry)
        path::String = strip(header[entry])
        if path == "NONE"
            return missing
        end
        name = find_file(data_directory, path, parent_name, exts)
        if !ismissing(name)
            return name
        end
    end
    missing
end

function find_file(dir, name, parent, extensions)
    if length(name) == 0
        return missing
    elseif match(r"%match%", name) !== nothing
        base = splitext(parent)[1]
        for ext in extensions
            testfile = joinpath(dir, base * ext)
            if isfile(testfile)
                return testfile
            end
        end
        @warn "Missing! Could not find file '%match%': tried $extensions"
        return missing
    elseif match(r"^none\b", name) !== nothing
        return missing
    end
    joinpath(dir, name)
end

end # module

using .OGIP
export OGIP, StandardOGIPConfig

function read_fits_header(path; hdu = 2)
    OGIP._read_fits_and_close(path) do f
        FITSIO.read_header(f[hdu])
    end
end

export read_fits_header
