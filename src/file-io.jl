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

function load_spectral_dataset(pha_path, rmf_path)
    fgrp = FITS(pha_path)
    raw_grp = DataFrame(fgrp[2])
    rmf = load_response_file(rmf_path)
    grp = augment_energy_channel(raw_grp, rmf)
    (grp, rmf)
end

export load_response_file, augment_energy_channel, load_spectral_dataset
