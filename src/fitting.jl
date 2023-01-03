export wrap_model, wrap_model_noflux

function wrap_model(
    model::AbstractSpectralModel,
    data::SpectralDataset{T};
    energy = energy_vector(data),
) where {T}
    fluxes = make_fluxes(energy, flux_count(model), T)
    frozen_params = get_value.(frozenparameters(model))
    ΔE = data.energy_bin_widths
    # pre-mask the response matrix to ensure channel out corresponds to the active data points
    R = fold_ancillary(data)[data.mask, :]
    # pre-allocate the output 
    outflux = zeros(T, length(ΔE))
    wrapped =
        (energy, params) -> begin
            invokemodel!(fluxes, energy, model, params, frozen_params)
            mul!(outflux, R, fluxes[1])
            @. outflux = outflux / ΔE
        end
    energy, wrapped
end

function wrap_model_noflux(
    model::AbstractSpectralModel,
    data::SpectralDataset{T};
    energy = energy_vector(data),
) where {T}
    ΔE = data.energy_bin_widths
    # pre-mask the response matrix to ensure channel out corresponds to the active data points
    R = fold_ancillary(data)[data.mask, :]
    # pre-allocate the output 
    wrapped = (energy, params) -> begin
        flux = invokemodel(energy, model, params)
        flux = (R * flux)
        @. flux = flux / ΔE
    end
    energy, wrapped
end
