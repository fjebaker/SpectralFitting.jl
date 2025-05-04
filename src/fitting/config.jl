struct ModelFittingCache{BackingCacheType}
    model_output::BackingCacheType
    objective_output::BackingCacheType
end

function ModelFittingCache(model_output::AbstractArray, objective_output::AbstractArray)
    ModelFittingCache(DiffCache(model_output), DiffCache(objective_output))
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(cache::ModelFittingCache))
    println(
        io,
        "ModelFittingCache[#model_ouput=$(size(cache.model_output.du)),#objective_output=$(size(cache.objective_output.du))]",
    )
end

function _make_model_cache(layout, model::AbstractSpectralModel, model_domain, objective)
    model_output = construct_objective_cache(layout, model, model_domain)
    objective_output = zeros(eltype(objective), (length(objective), 1))
    ModelFittingCache(model_output, objective_output)
end

function calculate_objective!(
    cache::ModelFittingCache,
    domain,
    model::AbstractSpectralModel{T},
    transformer!!,
) where {T<:Number}
    model_output = get_tmp(cache.model_output, T)
    objective_output = get_tmp(cache.objective_output, T)

    output = invokemodel!(model_output, domain, model)
    transformer!!(objective_output, domain, output)

    @views objective_output[:, 1]
end

struct DatasetFittingCache{V}
    model_domain::V
    objective_domain::V
    objective::V
    variance::V
    covariance::V
end

struct FittingConfig{
    T,
    Prob<:FittingProblem,
    Statistic,
    ParameterCacheType,
    BindingsType,
    ModelCacheType,
    DataCacheType,
    TransformerType,
}
    u0::Vector{FitParam{T}}
    prob::Prob
    stat::Statistic
    parameter_cache::ParameterCacheType
    parameter_bindings::BindingsType
    model_cache::ModelCacheType
    data_cache::DataCacheType
    transformers::TransformerType
end

fit_statistic(cfg::FittingConfig) = cfg.stat

function FittingConfig(prob::FittingProblem; stat = ChiSquared())
    # TODO: remove the bound parameters
    p_vector, bindings = parameter_vector_and_bindings(prob)

    free_mask = _make_free_mask(p_vector)
    v_vector = get_value.(p_vector)

    p_cache = ParameterCache(free_mask, DiffCache(v_vector), v_vector[.!free_mask])

    I = ((1:model_count(prob))...,)

    _res = map(I) do i
        model = prob.model.m[i]
        dataset = prob.data.d[i]

        layout =
            with_units(common_support(model, dataset), preferred_units(dataset, stat))

        model_domain = make_model_domain(layout, dataset)
        objective_domain = make_output_domain(layout, dataset)
        objective = make_objective(layout, dataset)
        objective_variance = make_objective_variance(layout, dataset)
        transformer!! = objective_transformer(layout, dataset)

        dataset_cache = DatasetFittingCache(
            model_domain,
            objective_domain,
            objective,
            objective_variance,
            inv.(objective_variance),
        )

        transformer!!,
        dataset_cache,
        _make_model_cache(layout, model, model_domain, objective)
    end

    all_transformers = map(i -> _res[i][1], I)
    all_data_cache = map(i -> _res[i][2], I)
    all_model_cache = map(i -> _res[i][3], I)

    FittingConfig(
        p_vector[free_mask],
        prob,
        stat,
        p_cache,
        bindings,
        all_model_cache,
        all_data_cache,
        all_transformers,
    )
end

function calculate_objective!(config::FittingConfig, parameters, i::Int)
    parameters = @views parameters[config.parameter_bindings[i]]
    model = remake_with_parameters(config.prob.model.m[i], parameters)
    calculate_objective!(
        config.model_cache[i],
        config.data_cache[i].model_domain,
        model,
        config.transformers[i],
    )
end

function measure_objective!(config::FittingConfig, u0, i::Int)
    all_parameters = update_free_parameters!(config.parameter_cache, u0)

    y = calculate_objective!(config, all_parameters, i)
    measure(
        fit_statistic(config),
        config.data_cache[i].objective,
        y,
        config.data_cache[i].variance,
    )
end

function calculate_objective!(config::FittingConfig, u0)
    all_parameters = update_free_parameters!(config.parameter_cache, u0)

    I = ((1:model_count(config.prob))...,)
    map(I) do i
        calculate_objective!(config, all_parameters, i)
    end
end

function _unpack_config(prob::FittingProblem; stat = ChiSquared(), kwargs...)
    config = FittingConfig(prob; stat = stat)
    kwargs, config
end

function supports_autodiff(cfg::FittingConfig)
    all(implementation(m) isa JuliaImplementation for m in cfg.prob.model.m)
end

function Base.show(io::IO, @nospecialize(config::FittingConfig))
    descr = "FittingConfig"
    print(io, descr)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(config::FittingConfig))
    descr = "FittingConfig"
    print(io, descr)
end

export FittingConfig, calculate_objective!
