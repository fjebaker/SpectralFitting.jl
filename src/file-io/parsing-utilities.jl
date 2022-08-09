function energy_to_channel(energy, channels, low_energy_bins, high_energy_bins)
    for (c, low_e, high_e) in zip(channels, low_energy_bins, high_energy_bins)
        if (energy ≥ low_e) && (high_e > energy)
            return convert(Int, c)
        end
    end
    # for type stability
    return -1
end

function build_matrix_response!(
    R,
    # from rm table
    rm_low_energy,
    rm_high_energy,
    F_chan,
    N_chan,
    matrix_lookup,
    # from energy table
    low_energy,
    high_energy,
    channels,
)
    for (i, (low_e, high_e, F, N)) in
        enumerate(zip(rm_low_energy, rm_high_energy, F_chan, N_chan))
        M = matrix_lookup[i]
        # find which channel
        av_energy = (high_e + low_e) / 2

        channel = energy_to_channel(av_energy, channels, low_energy, high_energy)
        if channel < 1
            @warn "No channel for response energy $(low_e) to $(high_e). Skipping."
            continue
        end

        index = 1
        for (first, len) in zip(F, N)
            if len == 0
                break
            end
            @views R[first:first+len-1, channel] .= M[index:index+len-1]
            index += len
        end
    end
end

function build_matrix_response(
    # from rm table
    rm_low_energy,
    rm_high_energy,
    F_chan,
    N_chan,
    matrix_lookup,
    # from energy table
    low_energy,
    high_energy,
    channels;
    T::Type = Float64,
)
    n = length(channels)
    R = spzeros(T, n, n)
    build_matrix_response!(
        R,
        # from rm table
        rm_low_energy,
        rm_high_energy,
        F_chan,
        N_chan,
        matrix_lookup,
        # from energy table
        low_energy,
        high_energy,
        channels,
    )
    R
end

function augmented_energy_channels(channels, rm::ResponseMatrix{T}) where {T}
    # just going to assume the channels line up
    N = length(channels)
    Emax = zeros(T, N)
    Emin = zeros(T, N)
    @inbounds for (i, c) in enumerate(channels)
        if c ≤ N
            index = findfirst(==(c), rm.channels)
            Emax[i] = rm.high_energy_bins[index]
            Emin[i] = rm.low_energy_bins[index]
        end
    end
    (Emin, Emax)
end
