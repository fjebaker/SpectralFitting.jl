
mutable struct SimulatedSpectrum{T,F} <: AbstractDataset
    input_domain::Vector{T}
    output_domain::Vector{T}
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

bin_widths(dataset::SimulatedSpectrum) = diff(dataset.output_domain)
plotting_domain(dataset::SimulatedSpectrum) =
    dataset.output_domain[1:end-1] .+ (bin_widths(dataset) ./ 2)
objective_units(dataset::SimulatedSpectrum) = dataset.units

function _printinfo(io::IO, spectrum::SimulatedSpectrum)
    dmin, dmax = prettyfloat.(extrema(spectrum.data))
    xmin, xmax = prettyfloat.(extrema(spectrum.domain))
    descr = """SimulatedSpectrum with $(length(spectrum.data)) channels:
      Units                 : $(spectrum.units)
      . Data (min/max)      : ($dmin, $dmax)
      . Domain (min/max)    : ($xmin, $xmax)
    """
    print(io, descr)
end

function simulate!(
    config::FittingConfig,
    p;
    simulate_distribution = Distributions.Poisson,
    rng = Random.default_rng(),
    exposure_time = 1e3,
    bin_width = nothing,
)
    config.objective .= _invoke_and_transform!(config.cache, config.domain, p)
    for (i, m) in enumerate(config.objective)
        distr = simulate_distribution(m * exposure_time)
        count = rand(rng, distr)
        config.objective[i] = count / (exposure_time)
        config.variance[i] = sqrt(abs(count)) / (exposure_time)
    end
end

function simulate!(
    output_domain::AbstractVector,
    conf::FittingConfig;
    seed = abs(rand(Int)),
    kwargs...,
)
    rng = Random.default_rng(seed)
    Random.seed!(rng, seed)
    simulate!(conf, get_value.(conf.parameters); rng = rng, kwargs...)
    SimulatedSpectrum(
        conf.domain,
        output_domain,
        conf.objective,
        # variance has already been sqrt'd
        conf.variance,
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
    input_domain = response_energy(response),
    output_domain = channel_energy(response),
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

    E = output_domain
    ΔE = diff(E)

    objective = zeros(eltype(E), length(E) - 1)
    variance = ones(eltype(objective), size(objective))

    domain = input_domain
    cache = SpectralCache(
        layout,
        model,
        domain,
        objective,
        _fold_transformer(T, layout, R, ΔE, E),
    )

    conf = FittingConfig(
        implementation(model),
        cache,
        _allocate_free_parameters(model),
        domain,
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
    input_domain = response_energy(response)
    output_domain = channel_energy(response)
    kw, conf = _make_simulation_fitting_config(
        model,
        response,
        ancillary;
        output_domain = output_domain,
        input_domain = input_domain,
        kwargs...,
    )
    simulate!(output_domain, conf; kw...)
end

function simulate(model::AbstractSpectralModel, dataset::AbstractDataset; kwargs...)
    kw, conf = _unpack_fitting_configuration(FittingProblem(model => dataset); kwargs...)
    simulate!(response_energy(dataset.data.response), conf; kw...)
end


export simulate
