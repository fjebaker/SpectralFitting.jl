@fastmath function integrate_over_flux!(callback, flux, energy)
    E1 = callback(first(energy))
    @inbounds for i in eachindex(flux)
        E2 = callback(energy[i+1])
        flux[i] = E2 - E1
        E1 = E2
    end
end
