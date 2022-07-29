struct ResponseMatrix{T,E}
    matrix::SparseMatrixCSC{T,Int64}
    ebins::E
end

foldflux(flux, rmf::ResponseMatrix) = rmf.matrix * flux
Base.:*(rmf::ResponseMatrix, flux) = rmf.matrix * flux

function energybins(rmf::ResponseMatrix{T}) where {T}
    energy = zeros(T, length(rmf.ebins.E_MAX) + 1)
    energy[1:end-1] .= rmf.ebins.E_MIN
    energy[end] = rmf.ebins.E_MAX[end]
    energy
end

function Base.show(
    io::IO,
    ::MIME{Symbol("text/plain")},
    rmf::ResponseMatrix{E,M},
) where {E,M}
    nchans = nrow(rmf.ebins)
    println(io, "ReponseMatrix with $nchans channels:")
    Base.print_array(io, rmf.matrix)
end

function energy_to_channel(energy, ebins)
    for j in eachrow(ebins)
        if (energy ≥ j.E_MIN) && (j.E_MAX > energy)
            return convert(Int64, j.CHANNEL)
        end
    end
    # for type stability
    return -1
end

function build_matrix_response!(R, table, ebins, matrix_lookup)
    for (i, row) in enumerate(eachrow(table))
        M = matrix_lookup[i]
        # find which channel
        av_energy = (row.ENERG_HI + row.ENERG_LO) / 2

        channel = energy_to_channel(av_energy, ebins)
        if channel < 1
            @warn "No channel for response energy $(row.ENERG_LO) to $(row.ENERG_HI). Skipping."
            continue
        end

        index = 1
        for (first, len) in zip(row.F_CHAN, row.N_CHAN)
            if len == 0
                break
            end
            @views R[first:first+len-1, channel] .= M[index:index+len-1]
            index += len
        end
    end
end

function build_matrix_response(table, ebins, matrix_lookup)
    n = nrow(ebins)
    R = spzeros(Float64, n, n)
    build_matrix_response!(R, table, ebins, matrix_lookup)
    R
end

function load_response_file(path)
    f = FITS(path)
    table = DataFrame(f[2])
    energy_bins = DataFrame(f[3])
    matrix_lookup = read(f[2], "MATRIX")

    R = build_matrix_response(table, energy_bins, matrix_lookup)
    ResponseMatrix(R, energy_bins)
end

# can we do this without copying?
function augment_energy_channel(grp, rmf::ResponseMatrix)
    joined = innerjoin(grp, rmf.ebins, on = :CHANNEL)
    joined.E_DELTA = joined.E_MAX .- joined.E_MIN
    joined
end

export ResponseMatrix, foldflux, load_response_file, augment_energy_channel, energybins
