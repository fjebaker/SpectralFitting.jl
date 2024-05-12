struct SimulatedSpectrum{T} <: AbstractDataset
    data::Vector{T}
    errors::Vector{T}
end

function simulate!(config::FittingConfig, p; simulate_distribution = Distributions.Normal)
    config.objective .= _invoke_and_transform!(config.cache, config.domain, p)
    for (i, m) in enumerate(config.objective)
        distr = simulate_distribution(m, sqrt(config.variance[i]))
        config.objective[i] = rand(distr)
    end
end

function simulate(prob::FittingProblem; kwargs...)
    kw, conf = _unpack_fitting_configuration(prob; kwargs...)
    simulate!(conf, conf.parameters; kw...)
    SimulatedSpectrum(conf.objective, sqrt.(conf.variance))
end

function simulate(model::AbstractSpectralModel, dataset::AbstractDataset; kwargs...)
    simulate(FittingProblem(model => dataset); kwargs...)
end


export simulate
