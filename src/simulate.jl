
mutable struct SimulatedSpectrum{T,F} <: AbstractDataset
    domain::Vector{T}
    data::Vector{T}
    errors::Vector{T}
    units::Union{Nothing,SpectralUnits.RateOrCount}
    transformer!!::F
    seed::Int
end

supports_contiguosly_binned(::Type{<:SimulatedSpectrum}) = true

function make_objective(::ContiguouslyBinned, dataset::SimulatedSpectrum)
    check_units_warning(dataset.units)
    dataset.data
end

function make_objective_variance(::ContiguouslyBinned, dataset::SimulatedSpectrum)
    check_units_warning(dataset.units)
    dataset.errors .^ 2
end

function make_model_domain(::ContiguouslyBinned, dataset::SimulatedSpectrum)
    dataset.domain
end

bin_widths(dataset::SimulatedSpectrum) = diff(dataset.domain)
plotting_domain(dataset::SimulatedSpectrum) = dataset.domain[1:end-1] .+ bin_widths(dataset)
objective_units(dataset::SimulatedSpectrum) = dataset.units

function _printinfo(io::IO, spectrum::SimulatedSpectrum)
    dmin, dmax = prettyfloat.(extrema(spectrum.data))
    descr = """SimulatedSpectrum:
      Units                 : $(spectrum.units)
      . Data (min/max)      : ($dmin, $dmax)
    """
    print(io, descr)
end

function simulate!(
    config::FittingConfig,
    p;
    simulate_distribution = Distributions.Normal,
    rng = Random.default_rng(),
)
    config.objective .= _invoke_and_transform!(config.cache, config.domain, p)
    for (i, m) in enumerate(config.objective)
        distr = simulate_distribution(m, sqrt(config.variance[i]))
        config.objective[i] = rand(rng, distr)
    end
end

function simulate!(conf::FittingConfig; seed = abs(rand(Int)), kwargs...)
    rng = Random.default_rng(seed)
    Random.seed!(rng, seed)
    simulate!(conf, get_value.(conf.parameters); rng = rng, kwargs...)
    SimulatedSpectrum(
        conf.domain,
        conf.objective,
        sqrt.(conf.variance),
        u"counts / (s * keV)",
        conf.cache.transformer!!,
        seed,
    )
end

function _make_simulation_fitting_config(
    model::AbstractSpectralModel,
    response::ResponseMatrix{T},
    ancillary;
    layout = ContiguouslyBinned(),
    kwargs...,
) where {T}
    if !supports(layout, model)
        throw("Model must support desired layout for simulation.")
    end

    R = if !isnothing(ancillary)
        fold_ancillary(response, ancillary)
    else
        response.matrix
    end

    E = response_energy(response)
    ΔE = diff(E)

    objective = zeros(eltype(E), length(E) - 1)
    variance = fill(1e-4, size(objective))
    cache =
        SpectralCache(layout, model, E, objective, _fold_transformer(T, layout, R, ΔE, E))

    conf = FittingConfig(
        implementation(model),
        cache,
        _allocate_free_parameters(model),
        E,
        objective,
        variance,
    )
    kwargs, conf
end

function simulate(
    model::AbstractSpectralModel,
    response::ResponseMatrix,
    ancillary::Union{Nothing,<:AncillaryResponse};
    kwargs...,
)
    kw, conf = _make_simulation_fitting_config(model, response, ancillary; kwargs...)
    simulate!(conf; kw...)
end

function simulate(model::AbstractSpectralModel, dataset::AbstractDataset; kwargs...)
    kw, conf = _unpack_fitting_configuration(FittingProblem(model => dataset); kwargs...)
    simulate!(conf; kw...)
end


export simulate
