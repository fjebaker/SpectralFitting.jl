export ResponseMatrix, fold_response, energy_bins

function ResponseMatrix(
    fits,
    mission::AbstractMission,
    ::Type{T},
)::ResponseMatrix{T} where {T}
    (matrix, rm_channels, E_lows, E_highs, rm_E_lows, rm_E_highs) = parse_rm_fits_file(mission, fits, T)
    ResponseMatrix(matrix, Int.(rm_channels), T.(E_lows), T.(E_highs), T.(rm_E_lows), T.(rm_E_highs))
end

# folding response

@fastmath fold_response(flux, rm::ResponseMatrix) = rm.matrix * flux

function fold_response(flux, energy, rm::ResponseMatrix)
    if length(rm.energy_bins_low) != length(energy) - 1
        out_flux = rebin_flux(flux, energy, rm)
        fold_response(out_flux, rm)
    else
        fold_response(flux, rm)
    end
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, rm::ResponseMatrix{T}) where {T}
    nchans = length(rm.channels)
    println(io, "ResponseMatrix with $nchans channels:")
    Base.print_array(io, rm.matrix)
end

function group_response(rm::ResponseMatrix{T}, grouping) where {T}

    return deepcopy(rm)

    # indices = grouping_to_indices(grouping)
    # N = length(indices) - 1
    # R = spzeros(T, (N, N))

    # energy_low = zeros(T, N)
    # energy_high = zeros(T, N)

    # grouping_indices_callback(indices) do (i, index1, index2)
    #     energy_low[i] = rm.energy_bins_low[index1]
    #     energy_high[i] = rm.energy_bins_high[index2]

    #     # iterate over row
    #     grouping_indices_callback(indices) do (j, index3, index4)
    #         R[i, j] = @views sum(rm.matrix[index1:index2, index3:index4])
    #         # R[i, j] = @views sum(rm.matrix[index3:index4, index1:index2])
    #     end

    #     # normalise row if non-zero
    #     Σr = sum(R[i, :])
    #     if Σr > 0.0
    #         @views R[i, :] .= R[i, :] ./ Σr
    #     end
    # end

    # # normalise

    # ResponseMatrix(R, collect(1:N), energy_low, energy_high)
end

function build_matrix_response!(
    R,
    # from rm table
    rmf_E_lows,
    rmf_E_highs,
    Fchan,
    Nchan,
    matrix_lookup,
    # from energy table
    low_energy,
    high_energy,
    channels
)
    for (i, (low_e, high_e, F, N)) in enumerate(zip(rmf_E_lows, rmf_E_highs, Fchan, Nchan))
        M = matrix_lookup[i]
        # find which channel
        av_energy = low_e

        channel = energy_to_channel(av_energy, channels, low_energy, high_energy)
        if channel < 0
            @warn "No channel for response energy $(low_e) to $(high_e). Skipping."
            continue
        end

        index = 1
        for (first, len) in zip(F, N)
            if len == 0
                break
            end
            @views R[first+1:first+len, i] .= M[index:index+len-1]
            index += len
        end
    end
end

# must return these four things
function parse_rm_fits_file(mission::AbstractMission, fits, T)
    # from matrix table
    rmf_E_lows = read(fits[2], "ENERG_LO")
    rmf_E_highs = read(fits[2], "ENERG_HI")
    Fchan = eachcol(read(fits[2], "F_CHAN"))
    Nchan = eachcol(read(fits[2], "N_CHAN"))
    matrix_lookup = read(fits[2], "MATRIX")
    # from energy table
    E_lows = read(fits[3], "E_MIN")
    E_highs = read(fits[3], "E_MAX")
    channels = read(fits[3], "CHANNEL") 

    N_rows = Int(read_header(fits[2])["DETCHANS"])
    N_cols = length(rmf_E_lows)
    matrix = spzeros(T, N_rows, N_cols)

    build_matrix_response!(
        matrix,
        # from rm table
        rmf_E_lows,
        rmf_E_highs,
        Fchan,
        Nchan,
        matrix_lookup,
        # from energy table
        E_lows,
        E_highs,
        channels
    )

    # normalise_rows!(matrix)

    return (matrix, channels, E_lows, E_highs, rmf_E_lows, rmf_E_highs)
end


function normalise_rows!(matrix)
    @views for i in 1:size(matrix, 2)
        ΣR = sum(matrix[i, :])
        if ΣR > 0.0
            matrix[i, :] .= matrix[i, :] / ΣR
        end
    end
end