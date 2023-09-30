export rebin_flux, construct_objective_cache, regroup!

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

construct_objective_cache(model::AbstractSpectralModel, domain::AbstractVector) = construct_objective_cache(preferred_support(model), model, domain)
construct_objective_cache(T::Type, model::AbstractSpectralModel, domain::AbstractVector) = construct_objective_cache(preferred_support(model), T, model, domain)
construct_objective_cache(layout::AbstractDataLayout, model::AbstractSpectralModel, domain::AbstractVector{T}) where {T} =
    construct_objective_cache(layout, T, model, domain)
function construct_objective_cache(layout::AbstractDataLayout, T::Type, model::AbstractSpectralModel, domain)
    N = if layout isa OneToOne
        length(domain)
    elseif layout isa ContiguouslyBinned
        length(domain) - 1
    else
        error("No domain construction for $(layout) defined.")
    end
    zeros(T, (N, objective_cache_count(model)))
end

function augmented_energy_channels(channels, other_channels, bins_low, bins_high)
    # TODO: just going to assume the channels line up
    N = length(channels)
    energies = zeros(eltype(bins_high), N + 1)
    for (i, c) in enumerate(channels)
        index = findnext(==(c), other_channels, i)
        if isnothing(index)
            error("Failed to calculate channel to energy mapping.")
        end
        if index > N
            break
        end
        if (i > 1) && !(energies[i] ≈ bins_low[index])
            @warn "Channel $index: misaligned $(energies[i-1]) != $(bins_low[index])! Data may not be contiguous."
        end
        energies[i] = bins_low[index]
        energies[i+1] = bins_high[index]
    end
    energies
end
