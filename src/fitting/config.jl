struct FittingConfig{ImplType,CacheType,StatT,ProbT,P,D,O}
    cache::CacheType
    stat::StatT
    prob::ProbT
    parameters::P
    model_domain::D
    output_domain::D
    objective::O
    variance::O
    covariance::O
    function FittingConfig(
        impl::AbstractSpectralModelImplementation,
        cache::C,
        stat::AbstractStatistic,
        prob::FP,
        params::P,
        model_domain::D,
        output_domain::D,
        objective::O,
        variance::O;
        covariance::O = inv.(variance),
    ) where {C<:AbstractFittingCache,FP,P,D,O}
        new{typeof(impl),C,typeof(stat),FP,P,D,O}(
            cache,
            stat,
            prob,
            params,
            model_domain,
            output_domain,
            objective,
            variance,
            covariance,
        )
    end
end

fit_statistic(::Type{<:FittingConfig{Impl,Cache,Stat}}) where {Impl,Cache,Stat} = Stat()
fit_statistic(::T) where {T<:FittingConfig} = fit_statistic(T)

function make_single_config(prob::FittingProblem, stat::AbstractStatistic)
    model = prob.model.m[1]
    dataset = prob.data.d[1]

    layout = with_units(common_support(model, dataset), preferred_units(dataset, stat))
    model_domain = make_model_domain(layout, dataset)
    output_domain = make_output_domain(layout, dataset)
    objective = make_objective(layout, dataset)
    variance = make_objective_variance(layout, dataset)
    params::Vector{paramtype(model)} = collect(filter(isfree, parameter_tuple(model)))
    cache = SpectralCache(
        layout,
        model,
        model_domain,
        objective,
        objective_transformer(layout, dataset),
    )
    FittingConfig(
        implementation(model),
        cache,
        stat,
        prob,
        params,
        model_domain,
        output_domain,
        objective,
        variance,
    )
end

function make_multi_config(prob::FittingProblem, stat::AbstractStatistic)
    impl =
        all(model -> implementation(model) isa JuliaImplementation, prob.model.m) ?
        JuliaImplementation() : XSPECImplementation()

    layout = common_support(prob.model.m..., prob.data.d...)

    variances = map(d -> make_objective_variance(layout, d), prob.data.d)
    # build index mappings for pulling out the data
    domains, domain_mapping = _build_domain_mapping(layout, prob.data)
    output_domains, output_domain_mapping = _build_output_domain_mapping(layout, prob.data)
    objectives, objective_mapping = _build_objective_mapping(layout, prob.data)
    parameters, parameter_mapping = _build_parameter_mapping(prob.model, prob.bindings)

    i::Int = 1
    caches = map(prob.model.m) do m
        c = SpectralCache(
            layout,
            m,
            domains[i],
            objectives[i],
            objective_transformer(layout, prob.data.d[i]),
            param_diff_cache_size = length(parameters),
        )
        i += 1
        c
    end

    all_objectives = reduce(vcat, objectives)

    cache = MultiModelCache(
        caches,
        DiffCache(similar(all_objectives)),
        domain_mapping,
        output_domain_mapping,
        objective_mapping,
        parameter_mapping,
    )
    FittingConfig(
        impl,
        cache,
        stat,
        prob,
        parameters,
        reduce(vcat, domains),
        reduce(vcat, output_domains),
        all_objectives,
        reduce(vcat, variances),
    )
end

function FittingConfig(prob::FittingProblem; stat = ChiSquared())
    config = if model_count(prob) == 1 && data_count(prob) == 1
        make_single_config(prob, stat)
    elseif model_count(prob) == data_count(prob)
        make_multi_config(prob, stat)
    elseif model_count(prob) < data_count(prob)
        error("Single model, many data not yet implemented.")
    else
        error("Multi model, single data not yet implemented.")
    end

    return config
end

function _unpack_config(prob::FittingProblem; stat = ChiSquared(), kwargs...)
    config = FittingConfig(prob; stat = stat)
    kwargs, config
end

function _f_objective(config::FittingConfig)
    function f!!(domain, parameters)
        _invoke_and_transform!(config.cache, domain, parameters)
    end
end

function paramtype(
    ::FittingConfig{ImplType,CacheType,StatT,ProbT,P},
) where {ImplType,CacheType,StatT,ProbT,P}
    T = eltype(P)
    K = if T <: FitParam
        paramtype(T)
    else
        T
    end
    Vector{K}
end
supports_autodiff(::FittingConfig{<:JuliaImplementation}) = true
supports_autodiff(::FittingConfig) = false

function Base.show(io::IO, @nospecialize(config::FittingConfig))
    descr = "FittingConfig"
    print(io, descr)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(config::FittingConfig))
    descr = "FittingConfig"
    print(io, descr)
end

export FittingConfig
