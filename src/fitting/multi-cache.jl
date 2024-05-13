struct MultiModelCache{K,N,CacheTypes<:Tuple,ParameterMappingType} <: AbstractFittingCache
    caches::CacheTypes
    all_outputs::K
    domain_mapping::NTuple{N,Int}
    output_domain_mapping::NTuple{N,Int}
    objective_mapping::NTuple{N,Int}
    parameter_mapping::ParameterMappingType
end

function _get_range(mapping::NTuple, i)
    m_start = i == 1 ? 1 : mapping[i-1] + 1
    m_end = mapping[i]
    (m_start, m_end)
end

function _invoke_and_transform!(cache::MultiModelCache, domain, params)
    all_outputs = get_tmp(cache.all_outputs, params)

    for (i, ch) in enumerate(cache.caches)
        p = @views params[cache.parameter_mapping[i]]

        domain_start, domain_end = _get_range(cache.domain_mapping, i)
        objective_start, objective_end = _get_range(cache.objective_mapping, i)

        d = @views domain[domain_start:domain_end]
        all_outputs[objective_start:objective_end] .= _invoke_and_transform!(ch, d, p)
    end

    all_outputs
end

function _build_parameter_mapping(model::FittableMultiModel, bindings)
    parameters = map(_allocate_free_parameters, model.m)
    parameters_counts = _accumulated_indices(map(length, parameters))

    all_parameters = reduce(vcat, parameters)

    parameter_mapping, remove = _construct_bound_mapping(bindings, parameters_counts)
    # remove duplicate parameters that are bound
    deleteat!(all_parameters, remove)

    all_parameters, parameter_mapping
end

function _build_mapping_length(f, itt::Tuple)
    values = map(f, itt)
    mapping = _accumulated_indices(map(length, values))
    values, mapping
end

_build_objective_mapping(layout::AbstractDataLayout, dataset::FittableMultiDataset) =
    _build_mapping_length(i -> make_objective(layout, i), dataset.d)
_build_domain_mapping(layout::AbstractDataLayout, dataset::FittableMultiDataset) =
    _build_mapping_length(i -> make_model_domain(layout, i), dataset.d)
_build_output_domain_mapping(layout::AbstractDataLayout, dataset::FittableMultiDataset) =
    _build_mapping_length(i -> make_output_domain(layout, i), dataset.d)

function FittingConfig(prob::FittingProblem)
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
        parameters,
        reduce(vcat, domains),
        reduce(vcat, output_domains),
        all_objectives,
        reduce(vcat, variances),
    )
end

function finalize(
    config::FittingConfig{Impl,<:MultiModelCache},
    params;
    statistic = ChiSquared(),
    σparams = nothing,
) where {Impl}
    domain = config.model_domain
    cache = config.cache
    results = map(enumerate(cache.caches)) do (i, ch)
        p = @views params[cache.parameter_mapping[i]]
        σp = @views isnothing(σparams) ? nothing : σparams[cache.parameter_mapping[i]]

        domain_start, domain_end = _get_range(cache.domain_mapping, i)
        objective_start, objective_end = _get_range(cache.objective_mapping, i)

        d = @views domain[domain_start:domain_end]

        output = _invoke_and_transform!(ch, d, p)

        chi2 = measure(
            statistic,
            config.objective[objective_start:objective_end],
            output,
            config.variance[objective_start:objective_end],
        )
        (; chi2, p, σp)
    end

    unc = getindex.(results, :σp)
    unc_or_nothing = if any(isnothing, unc)
        nothing
    else
        unc
    end
    MultiFittingResult(
        getindex.(results, :chi2),
        getindex.(results, :p),
        unc_or_nothing,
        config,
    )
end
