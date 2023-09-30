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
    output = zeros(eltype(input), length(target_bins_high))
    downsample_rebin!(output, input, current_bins, target_bins_high)
end


"""
```
             [--i------] 
|   |   |   |   |   |   |   |   |   |  current
   |         |         |         |     target
             j
```
"""
function _downsample_rebin!(output, target, input, current)
    @assert length(output) == length(target) - 1
    @assert length(input) == length(current) - 1

    # ensure everything is zeroed
    @. output = 0

    i_start = searchsortedfirst(current, target[1])

    ω_start = _straddle_ratio(target[1], current, i_start)
    output[1] += (1 - ω_start) * input[i_start - 1]

    # track slices
    start::Int = i_start
    stop::Int = i_start
    
    # loop over each output bin
    # h is the high of the bin in target
    for (j, h) in @views enumerate(target[2:end])
        # find the first
        i = @views searchsortedfirst(current[start:end], h) + (start - 1)
        if i > lastindex(current) + 1
            break # edge case
        end
        stop = i

        ω = _straddle_ratio(h, current, stop)
        output[j] += @views sum(input[start:stop-2]) + ω * input[stop-1]
        # carry over 
        if j < lastindex(output)
            output[j+1] += (1 - ω) * input[stop-1]
        end

        start = stop
    end
end

"""
```
                      Δ
|              |---------------|    
 |        |          | h
               [-----] τ
```
"""
function _straddle_ratio(h, current, stop)
    Δ = current[stop] - current[stop - 1] # current width
    τ = h - current[stop - 1]
    τ / Δ
end

function downsample_rebin!(output, input, current_bins, target_bins)
    # ensure we are down-sampling, and not up-sampling
    if (length(current_bins) ≤ length(target_bins))
        throw(
            "Rebinning must down-sample to fewer destination bins than source bins. Use interpolation methods for up-sampling.",
        )
    end

    @assert length(output) == length(target_bins) - 1
    @assert length(input) == length(current_bins) - 1

    N = length(target_bins)

    # handle lower edge
    i_start = findfirst(>(first(target_bins)), current_bins)
    view_current_bins = @view current_bins[i_start:end]

    # find first dest energy that is greater than first src energy
    i = findfirst(>(first(view_current_bins)), target_bins)
    # view into target bins of interest
    trunc_bin_high = @view(target_bins[i:end])

    start = stop = 1
    for (fi, Ehigh) in zip(i:N, trunc_bin_high)
        # find lower and upper limit index
        stop = findnext(>(Ehigh), view_current_bins, stop)
        # break if no energy is higher
        isnothing(stop) && break

        # sum the ones that are between immediately
        output[fi] += @views sum(input[start:stop-2])

        # deal with edge bin by calculating overlap
        ΔE = view_current_bins[stop] - view_current_bins[stop-1]
        # amount of the bin in the current dest bin
        δE = Ehigh - view_current_bins[stop-1]
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

construct_objective_cache(model::AbstractSpectralModel, domain::AbstractVector) =
    construct_objective_cache(preferred_support(model), model, domain)
construct_objective_cache(T::Type, model::AbstractSpectralModel, domain::AbstractVector) =
    construct_objective_cache(preferred_support(model), T, model, domain)
construct_objective_cache(
    layout::AbstractDataLayout,
    model::AbstractSpectralModel,
    domain::AbstractVector{T},
) where {T} = construct_objective_cache(layout, T, model, domain)
function construct_objective_cache(
    layout::AbstractDataLayout,
    T::Type,
    model::AbstractSpectralModel,
    domain,
)
    construct_objective_cache(layout, T, length(domain), objective_cache_count(model))
end
function construct_objective_cache(layout::AbstractDataLayout, T::Type, N::Int, size::Int)
    n = if layout isa OneToOne
        N
    elseif layout isa ContiguouslyBinned
        N - 1
    else
        error("No domain construction for $(layout) defined.")
    end
    zeros(T, (n, size))
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
