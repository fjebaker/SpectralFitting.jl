export rebin_flux, make_flux, make_fluxes

rebin_flux(flux, curr_energy, dest_energy_bins::AbstractVector) =
    rebin_flux(flux, curr_energy, first(dest_energy_bins), @view(dest_energy_bins[2:end]))
rebin_flux(flux, curr_energy, rm::ResponseMatrix) =
    rebin_flux(flux, curr_energy, first(rm.energy_bins_low), rm.energy_bins_high)

function rebin_flux(flux, curr_energy, E_min::Number, energy_bins_high::AbstractVector)
    # ensure we are sub-sampling, and not up-sampling
    @assert (lenngth(curr_energy) ≤ length(energy_bins_high) + 1)

    N = length(energy_bins_high)
    out_flux = zeros(eltype(flux), N)

    # find the initial index into energy
    current_i = findfirst(>(E_min), curr_energy)

    for (flux_index, E_max_bins) in enumerate(energy_bins_high)
        # find where energy is outside of bin
        next_i = findnext(>(E_max_bins), curr_energy, current_i+1)
        isnothing(next_i) && break

        # sum everything that is definitely inbetween
        @views out_flux[flux_index] += sum(flux[current_i:next_i-2])

        # calculate overlap flux: width of energy bin
        ΔE = curr_energy[next_i] - curr_energy[next_i-1]
        # width of bin in current flux bin
        δE = E_max_bins - curr_energy[next_i-1]
        ratio = (δE / ΔE)
        @views out_flux[flux_index] += ratio * flux[next_i-1]

        # if not at end of flux array, carry over
        if flux_index < N
            out_flux[flux_index+1] += (1 - ratio) * flux[next_i-1]
        end
        current_i = next_i
    end

    out_flux
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
    for i in 1:length(indices)-1
        index1 = indices[i]
        index2 = indices[i+1] - 1
        func((i, index1, index2))
    end
    indices
end