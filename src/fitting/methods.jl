export AbstractFittingAlgorithm, LevenbergMarquadt, fit, fit!

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

function configuration(prob::FittingProblem; kwargs...)
    kw, config = _unpack_config(prob; kwargs...)
    if length(kw) > 0
        throw("Unknown keyword arguments: $(kw)")
    end
    config
end

function fit(prob::FittingProblem, args...; kwargs...)
    method_kwargs, config = _unpack_config(prob; kwargs...)
    @inline fit(config, args...; method_kwargs...)
end

function fit(
    config::FittingConfig,
    alg::LevenbergMarquadt;
    verbose = false,
    max_iter = 1000,
    method_kwargs...,
)
    @assert fit_statistic(config) == ChiSquared() "Least squares only for χ2 statistics."

    lsq_result = _lsq_fit(
        _f_objective(config),
        config.model_domain,
        config.objective,
        config.covariance,
        config.parameters,
        alg;
        verbose = verbose,
        max_iter = max_iter,
        autodiff = supports_autodiff(config) ? :forward : :finite,
        method_kwargs...,
    )
    params = LsqFit.coef(lsq_result)
    σ = try
        LsqFit.standard_errors(lsq_result)
    catch e
        @warn "No parameter uncertainty estimation due to error: $e"
        nothing
    end

    y = _f_objective(config)(config.model_domain, params)
    chi2 = measure(fit_statistic(config), config.objective, y, config.variance)
    finalize(config, params, chi2; σparams = σ)
end

function fit(
    config::FittingConfig,
    optim_alg;
    verbose = false,
    autodiff = _determine_ad_backend(config),
    method_kwargs...,
)
    objective = _f_wrap_objective(fit_statistic(config), config)
    u0 = get_value.(config.parameters)

    if !(autodiff isa Optimization.SciMLBase.NoAD) && (!supports_autodiff(config))
        error("Model does not support automatic differentiation.")
    end

    lb, ub = _determine_bounds(config, autodiff)

    # build problem and solve
    opt_f = Optimization.OptimizationFunction(objective, autodiff)
    # todo: something is broken with passing the boundaries
    opt_prob =
        Optimization.OptimizationProblem(opt_f, u0, config.model_domain; lb = lb, ub = ub)

    sol = Optimization.solve(opt_prob, optim_alg; method_kwargs...)

    # TODO: temporary fix for type instabilities in Optimizations.jl
    new_pars::paramtype(config) = sol.u
    final_stat = objective(new_pars, config.model_domain)
    finalize(config, new_pars, final_stat)
end

function fit!(prob::FittingProblem, args...; kwargs...)
    result = fit(prob, args...; kwargs...)
    update_model!(prob.model, result)
    result
end

function _determine_bounds(config, ::A) where {A}
    if A <: Optimization.SciMLBase.NoAD
        nothing, nothing
    else
        get_lowerlimit.(config.parameters), get_upperlimit.(config.parameters)
    end
end

function _determine_ad_backend(config)
    if supports_autodiff(config)
        Optimization.AutoForwardDiff()
    else
        Optimization.SciMLBase.NoAD()
    end
end
