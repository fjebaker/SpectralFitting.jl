export rebin_flux, make_flux, make_fluxes

function rebin_flux(flux, current_energy, dest_energy_bins::AbstractVector)
    downsample_rebin(
        flux,
        current_energy,
        #@view(dest_energy_bins[1:end-1]),
        @view(dest_energy_bins[2:end])
    )
end

function rebin_flux(flux, current_energy, rm::ResponseMatrix)
    downsample_rebin(flux, current_energy, rm.energy_bins_high)
end

function downsample_rebin(input, current_bins, target_bins_high)
    # ensure we are down-sampling, and not up-sampling
    if (length(current_bins) ≤ length(target_bins_high) + 1)
        throw(
            "Rebinning must down-sample to fewer destination bins than source bins. Use interpolation methods for up-sampling.",
        )
    end

    N = length(target_bins_high)
    output = zeros(eltype(input), N)

    # find first dest energy that is greater than first src energy
    i = findfirst(>(first(current_bins)), target_bins_high)
    # view into target bins of interest
    trunc_bin_high = @view(target_bins_high[i:end])

    start = stop = 1
    for (fi, Ehigh) in zip(i:N, trunc_bin_high)
        # find lower and upper limit index
        stop = findnext(>(Ehigh), current_bins, stop)
        # break if no energy is higher
        isnothing(stop) && break

        # sum the ones that are between immediately
        output[fi] += @views sum(input[start:stop-2])

        # deal with edge bin by calculating overlap
        ΔE = current_bins[stop] - current_bins[stop-1]
        # amount of the bin in the current dest bin
        δE = Ehigh - current_bins[stop-1]
        ratio = (δE / ΔE)
        output[fi] += ratio * input[stop-1]

        # if not at end of input array, carry over the rest
        if fi < N
            output[fi+1] += (1 - ratio) * input[stop-1]
        end
        start = stop
    end
    output
end

function get_energy_bins(x, T::Type)
    energy = zeros(T, length(x.energy_bins_low) + 1)
    energy[1:end-1] .= x.energy_bins_low
    energy[end] = x.energy_bins_high[end]
    energy
end

make_flux(energy::AbstractVector{T}) where {T} = make_flux(T, length(energy) - 1)
make_flux(T::Type, length::Number) = zeros(T, length)

function make_fluxes(energy, N::Int)
    make_fluxes(energy, N, eltype(energy))
end

function make_fluxes(energy, N::Int, T::Type)
    flux = make_flux(T, length(energy) - 1)
    fluxes = typeof(flux)[flux]
    for _ = 1:N-1
        push!(fluxes, deepcopy(flux))
    end
    fluxes
end

function make_dual_fluxes(energy, N::Int)
    flux = make_flux(eltype(energy), length(energy) - 1)
    d_flux = dualcache(flux)
    d_fluxes = typeof(d_flux)[d_flux]
    for _ = 1:N-1
        push!(d_fluxes, deepcopy(d_flux))
    end
    d_fluxes
end

function augmented_energy_channels(channels, rm::ResponseMatrix{T}) where {T}
    # just going to assume the channels line up
    N = length(channels)
    Emax = zeros(T, N)
    Emin = zeros(T, N)
    @inbounds for (i, c) in enumerate(channels)
        if c ≤ N
            index = findfirst(==(c), rm.channels)
            Emax[i] = rm.energy_bins_high[index]
            Emin[i] = rm.energy_bins_low[index]
        end
    end
    (Emin, Emax)
end

function energy_to_channel(E, channels, energy_bins_low, energy_bins_high, start)
    @views energy_to_channel(
        E,
        channels[start:end],
        energy_bins_low[start:end],
        energy_bins_high[start:end],
    )
end
function energy_to_channel(E, channels, energy_bins_low, energy_bins_high)
    for (c, low_e, high_e) in zip(channels, energy_bins_low, energy_bins_high)
        if (E ≥ low_e) && (high_e > E)
            return convert(Int, c)
        end
    end
    # for type stability
    return -1
end

function grouping_to_indices(grouping)
    indices = Int[]
    for (i, g) in enumerate(grouping)
        if g == 1
            push!(indices, i)
        end
    end
    push!(indices, length(grouping))
    indices
end

function grouping_indices_callback(func, indices)
    for i = 1:length(indices)-1
        index1 = indices[i]
        index2 = indices[i+1] - 1
        func((i, index1, index2))
    end
end
