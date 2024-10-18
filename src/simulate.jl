
mutable struct SimulatedSpectrum{T,F} <: AbstractDataset
    model_domain::Vector{T}
    output_domain::Vector{T}
    data::Vector{T}
    variance::Vector{T}
    units::Union{Nothing,SpectralUnits.RateOrCount}
    transformer!!::F
    seed::Int
end

supports(::Type{<:SimulatedSpectrum}) = (ContiguouslyBinned(),)

function make_objective(::ContiguouslyBinned, dataset::SimulatedSpectrum)
    dataset.data
end

function make_objective_variance(::ContiguouslyBinned, dataset::SimulatedSpectrum)
    dataset.variance
end

make_model_domain(::ContiguouslyBinned, dataset::SimulatedSpectrum) = dataset.model_domain
make_output_domain(::ContiguouslyBinned, dataset::SimulatedSpectrum) = dataset.output_domain

objective_transformer(::ContiguouslyBinned, dataset::SimulatedSpectrum) =
    dataset.transformer!!

bin_widths(dataset::SimulatedSpectrum) = diff(dataset.output_domain)
plotting_domain(dataset::SimulatedSpectrum) =
    dataset.output_domain[1:(end-1)] .+ (bin_widths(dataset) ./ 2)
objective_units(dataset::SimulatedSpectrum) = dataset.units

function _printinfo(io::IO, spectrum::SimulatedSpectrum)
    dmin, dmax = prettyfloat.(extrema(spectrum.data))
    xmin, xmax = prettyfloat.(extrema(spectrum.model_domain))
    omin, omax = prettyfloat.(extrema(spectrum.output_domain))
    descr = """SimulatedSpectrum with $(length(spectrum.data)) channels:
      Units                 : $(spectrum.units)
      . Data (min/max)      : ($dmin, $dmax)
      . Domain (min/max)    : ($xmin, $xmax)
      . Out Dmn. (min/max)  : ($omin, $omax)
    """
    print(io, descr)
end

function simulate!(
    config::FittingConfig,
    p;
    simulate_distribution = Distributions.Poisson,
    rng = Random.default_rng(),
    exposure_time = 1e5,
    counting_variance = true,
)
    config.objective .= _invoke_and_transform!(config.cache, config.model_domain, p)
    for (i, m) in enumerate(config.objective)
        distr = simulate_distribution(m * exposure_time)
        count = rand(rng, distr)
        config.objective[i] = count / exposure_time

        # how to propagate the variance
        if counting_variance
            config.variance[i] = count / exposure_time^2
        end
    end
end

function simulate!(conf::FittingConfig; seed = abs(rand(Int)), kwargs...)
    rng = Random.default_rng(seed)
    Random.seed!(rng, seed)
    simulate!(conf, get_value.(conf.parameters); rng = rng, kwargs...)
    SimulatedSpectrum(
        conf.model_domain,
        conf.output_domain,
        conf.objective,
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
    output_domain = folded_energy(response),
    stat = ChiSquared(),
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

    ΔE = diff(output_domain)

    objective = zeros(eltype(output_domain), length(output_domain) - 1)
    variance = ones(eltype(objective), size(objective))

    cache = SpectralCache(
        layout,
        model,
        input_domain,
        objective,
        _fold_transformer(T, one(eltype(ΔE)), layout, R, ΔE, input_domain),
    )

    free_params = collect(filter(isfree, parameter_tuple(model)))

    conf = FittingConfig(
        implementation(model),
        cache,
        stat,
        nothing,
        free_params,
        input_domain,
        output_domain,
        objective,
        variance,
    )
    kwargs, conf
end

function simulate(
    model::AbstractSpectralModel,
    response::ResponseMatrix,
    ancillary::Union{Nothing,<:AncillaryResponse} = nothing;
    kwargs...,
)
    kw, conf = _make_simulation_fitting_config(model, response, ancillary; kwargs...)
    simulate!(conf; kw...)
end

function simulate(model::AbstractSpectralModel, dataset::AbstractDataset; kwargs...)
    kw, conf = _unpack_config(FittingProblem(model => dataset); kwargs...)
    simulate!(conf; kw...)
end

"""
    simulate(model::AbstractSpectralModel; kwargs...)

Returns an [`InjectiveDataset`](@ref) with a realisation of the model.
Can also add noise to the objective, but not to the domain.

The `kwargs` are:
- `seed`: if not `nothing` used to set the PRNG seed.
- `var`: variance of the data (assuming mean of 0)
"""
function simulate(
    domain::AbstractVector,
    model::AbstractSpectralModel;
    seed::Union{Nothing,Int} = nothing,
    simulate_distribution = Distributions.Normal,
    var = 0.1,
)
    _seed::Int = isnothing(seed) ? time_ns() : seed
    rng = Random.default_rng()
    Random.seed!(rng, _seed)

    flux = invokemodel(domain, model)

    realisation = map(flux) do f
        rand(rng, simulate_distribution(f, sqrt(var)))
    end
    variances = fill(var, size(realisation))
    BinnedData(domain, realisation; codomain_variance = variances)
end


export simulate
