export wrap_model

function wrap_model(
    model::AbstractSpectralModel,
    data::SpectralDataset{T};
    energy = energy_vector(data),
) where {T}
    fluxes = make_fluxes(energy, flux_count(model), T)
    frozen_params = get_value.(get_frozen_model_params(model))
    ΔE = data.energy_bin_widths
    # pre-mask the response matrix to ensure channel out corresponds to the active data points
    R = data.response.matrix[data.mask, :]
    # pre-allocate the output 
    outflux = zeros(T, length(ΔE))
    (energy, params) -> begin
        invokemodel!(fluxes, energy, model, params, frozen_params)
        mul!(outflux, R, fluxes[1])
        @. outflux = outflux / ΔE
    end
end
