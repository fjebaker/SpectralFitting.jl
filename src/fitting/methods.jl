export wrap_model, wrap_model_noflux
export AbstractFittingAlgorithm, LevenbergMarquadt, fit

abstract type AbstractFittingAlgorithm end

struct LevenbergMarquadt{T} <: AbstractFittingAlgorithm
    λ_inc::T
    λ_dec::T
end
LevenbergMarquadt(; λ_inc = 10.0, λ_dec = 0.1) = LevenbergMarquadt(λ_inc, λ_dec)

function _lsq_fit(
    f,
    x,
    y,
    cov,
    parameters,
    alg;
    verbose = false,
    max_iter = 1000,
    kwargs...,
)
    LsqFit.curve_fit(
        f,
        x,
        y,
        cov,
        get_value.(parameters);
        lower = get_lowerlimit.(parameters),
        upper = get_upperlimit.(parameters),
        lambda_increase = alg.λ_inc,
        lambda_decrease = alg.λ_dec,
        show_trace = verbose,
        maxIter = max_iter,
        kwargs...,
    )
end

struct FittingConfig{U,X,Y,V}
    u::U
    x::X
    y::Y
    variance::V
    covariance::V
    autodiff::Bool
end

FittingConfig(u, x, y, variance; autodiff = false, covariance = inv.(variance)) =
    FittingConfig(u, x, y, variance, covariance, autodiff)

function _unpack_fitting_configuration(prob)
    if model_count(prob) == 1 && data_count(prob) == 1
        # get first/only model and data
        model = prob.model.m[1]
        data = prob.data.d[1]

        # wrap the model
        f = _lazy_folded_invokemodel(model, data)
        u = freeparameters(model)
        x = domain_vector(data)
        y = target_vector(data)
        variance = target_variance(data)
        autodiff = implementation(model) isa JuliaImplementation

        bundler_single(res) = bundle_result(res, model, f, x, y, variance)
        config = FittingConfig(u, x, y, variance; autodiff = autodiff)
        return f, config, bundler_single
    elseif model_count(prob) == data_count(prob)
        f, u, state = assemble_multimodel(prob)
        x = domain_vector(prob.data)
        y = target_vector(prob.data)
        variance = target_variance(prob.data)

        bundler_multi(res) = bundle_multiresult(res, prob.model, x, y, variance, state)
        config = FittingConfig(u, x, y, variance; autodiff = state.autodiff)
        return f, config, bundler_multi
    elseif model_count(prob) < data_count(prob)
        error("Single model, many data not yet implemented.")
    else
        error("Multi model, single data not yet implemented.")
    end
end

function fit(
    prob::FittingProblem,
    alg::LevenbergMarquadt;
    verbose = false,
    max_iter = 1000,
    kwargs...,
)
    f, config, bundler = _unpack_fitting_configuration(prob)
    lsq_result = _lsq_fit(
        f,
        config.x,
        config.y,
        config.covariance,
        config.u,
        alg;
        verbose = verbose,
        max_iter = max_iter,
        autodiff = config.autodiff ? :forward : :finite,
        kwargs...,
    )
    bundler(LsqFit.coef(lsq_result))
end

function wrap_model(
    model::AbstractSpectralModel,
    data::SpectralDataset{T};
    energy = domain_vector(data),
) where {T}
    fluxes = make_fluxes(model, energy)
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
    energy = domain_vector(data),
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

function wrap_Optimization(
    model::AbstractSpectralModel,
    data::SpectralDataset{T};
    energy = domain_vector(data),
    target = data.rate,
    variance = data.rateerror .^ 2,
) where {T}
    ΔE = data.energy_bin_widths
    # pre-mask the response matrix to ensure channel out corresponds to the active data points
    R = fold_ancillary(data)[data.mask, :]
    n = length(target)
    # pre-allocate the output 
    wrapped = (params, energy) -> begin
        flux = invokemodel(energy, model, params)
        flux = (R * flux)
        @. flux = flux / ΔE
        χ2_from_ŷyvar(flux, target, variance)
        # l = @. flux - target + target * ( log(target) - log(flux) )
        # -l
        # l, flux
    end
    energy, wrapped
end

χ2_from_ŷyvar(ŷ, y, variance) = sum(@.((y - ŷ)^2 / variance))

function χ2(model::Function, params, data::SpectralDataset; energy = domain_vector(data))
    ŷ = model(energy, params)
    χ2_from_ŷyvar(ŷ, data.rate, data.rateerror .^ 2)
end
