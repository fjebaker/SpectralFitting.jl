export rebin_flux, make_flux, make_fluxes, domain_vector, regroup!

function rebin_flux(flux, current_energy, dest_energy_bins::AbstractVector)
    downsample_rebin(
        flux,
        current_energy,
        #@view(dest_energy_bins[1:end-1]),
        @view(dest_energy_bins[2:end])
    )
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

make_flux(m::AbstractSpectralModel, e::AbstractVector) = make_flux(eltype(e), m, e)
make_flux(T::Type, m::AbstractSpectralModel, e::AbstractVector) =
    make_flux(T, length(e) + Δoutput_length(m))
make_flux(T::Type, n::Int) = zeros(T, n)

make_fluxes(model::AbstractSpectralModel, domain::AbstractVector{T}) where {T} =
    make_fluxes(T, model, domain)
function make_fluxes(T, model::AbstractSpectralModel, domain)
    zeros(T, (length(domain) - 1, flux_count(model)))
end

function augmented_energy_channels(channels, other_channels, bins_high, bins_low)
    # TODO: just going to assume the channels line up
    N = length(channels)
    Emax = zeros(eltype(bins_high), N)
    Emin = zeros(eltype(bins_high), N)
    @inbounds for (i, c) in enumerate(channels)
        if c ≤ N
            index = findfirst(==(c), other_channels)
            Emax[i] = bins_high[index]
            Emin[i] = bins_low[index]
        end
    end
    (Emin, Emax)
end
