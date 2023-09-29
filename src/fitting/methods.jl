export AbstractFittingAlgorithm, LevenbergMarquadt, fit

abstract type AbstractFittingAlgorithm end
abstract type AbstractStatistic end

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

function _unpack_multi_model_multi_data(prob; subtract_background = true)
    f, u, state = assemble_multimodel(prob, subtract_background = subtract_background)
    x = domain_vector(prob.data)
    y = background_subtracted_target_vector(prob.data, subtract_background)
    variance = background_subtracted_target_variance(prob.data, subtract_background)

    bundler(res) = bundle_multiresult(res, prob.model, x, y, variance, state)
end

function _unpack_fitting_configuration(prob; kwargs...)
    config = if model_count(prob) == 1 && data_count(prob) == 1
        FittingConfig(prob.model.m[1], prob.data.d[1])
    elseif model_count(prob) == data_count(prob)
        FittingConfig(prob)
    elseif model_count(prob) < data_count(prob)
        error("Single model, many data not yet implemented.")
    else
        error("Multi model, single data not yet implemented.")
    end

    kwargs, config
end

function fit(
    prob::FittingProblem,
    alg::LevenbergMarquadt;
    verbose = false,
    max_iter = 1000,
    kwargs...,
)
    method_kwargs, config = _unpack_fitting_configuration(prob; kwargs...)
    lsq_result = _lsq_fit(
        _f_objective(config),
        config.domain,
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
    finalize(config, params)
end

function fit(
    prob::FittingProblem,
    stat::AbstractStatistic,
    optim_alg;
    verbose = false,
    autodiff = nothing,
    kwargs...,
)
    method_kwargs, f, config = _unpack_fitting_configuration(prob; kwargs...)
    objective = wrap_objective(stat, f, config)
    u0 = get_value.(config.u)
    lower = get_lowerlimit.(config.u)
    upper = get_upperlimit.(config.u)

    # determine autodiff
    if !((isnothing(autodiff)) || (autodiff isa Optimization.SciMLBase.NoAD)) &&
       !supports_autodiff(config)
        error("Model does not support automatic differentiation.")
    end
    _autodiff = if supports_autodiff(config) && isnothing(autodiff)
        Optimization.AutoForwardDiff()
    elseif !isnothing(autodiff)
        autodiff
    else
        Optimization.SciMLBase.NoAD()
    end

    # build problem and solve
    opt_f = Optimization.OptimizationFunction{false}(objective, _autodiff)
    # todo: something is broken with passing the boundaries
    opt_prob = Optimization.OptimizationProblem{false}(opt_f, u0, config.x)
    sol = Optimization.solve(opt_prob, optim_alg; method_kwargs...)
    sol.u
end
