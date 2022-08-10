export rebin_flux, make_flux, make_fluxes

rebin_flux(flux, curr_energy, dest_energy_bins::AbstractVector) =
    rebin_flux(flux, curr_energy, first(dest_energy_bins), @view(dest_energy_bins[2:end]))
rebin_flux(flux, curr_energy, rm::ResponseMatrix) =
    rebin_flux(flux, curr_energy, first(rm.low_energy_bins), rm.high_energy_bins)

function rebin_flux(flux, curr_energy, E_min::Number, high_energy_bins::AbstractVector)
    N = length(high_energy_bins)
    out_flux = zeros(eltype(flux), N)

    # find the initial index into energy
    current_i = findfirst(>(E_min), curr_energy)

    for (flux_index, E_max_bins) in enumerate(high_energy_bins)
        # find where energy is outside of bin
        next_i = findnext(>(E_max_bins), curr_energy, current_i + 1)
        isnothing(next_i) && break

        # sum everything inbetween
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
