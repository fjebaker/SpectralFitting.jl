
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

function simulate(prob::FittingProblem; seed = abs(randint()), kwargs...)
    kw, conf = _unpack_fitting_configuration(prob; kwargs...)
    rng = Random.default_rng(seed)
    Random.seed!(rng, seed)
    simulate!(conf, get_value.(conf.parameters); rng = rng, kw...)
    SimulatedSpectrum(
        conf.domain,
        conf.objective,
        sqrt.(conf.variance),
        nothing,
        conf.cache.transformer!!,
        seed,
    )
end

function simulate(model::AbstractSpectralModel, dataset::AbstractDataset; kwargs...)
    simulate(FittingProblem(model => dataset); kwargs...)
end


export simulate
