measure(stat::AbstractStatistic, y, ŷ, σ²) = error("Implement missing for $(typeof(stat))")

struct ChiSquared <: AbstractStatistic end
measure(::ChiSquared, y, ŷ, σ²) = sum(@.((y - ŷ)^2 / σ²))

struct Cash <: AbstractStatistic end

function _f_wrap_objective(stat::AbstractStatistic, config::FittingConfig)
    f = _f_objective(config)
    function _objective(parameters, domain)
        ŷ = f(domain, parameters)
        measure(stat, config.objective, ŷ, config.variance)
    end
end

export ChiSquared, Cash
