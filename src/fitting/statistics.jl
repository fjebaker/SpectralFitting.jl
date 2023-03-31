measure(stat::AbstractStatistic, y, ŷ, σ²) = error("Implement missing for $(typeof(stat))")

struct ChiSquared <: AbstractStatistic end
measure(::ChiSquared, y, ŷ, σ²) = sum(@.((y - ŷ)^2 / σ²))

struct Cash <: AbstractStatistic end

function wrap_objective(stat::AbstractStatistic, config::FittingConfig)
    function _objective(u, x)
        ŷ = config.f(x, u)
        measure(stat, config.y, ŷ, config.variance)
    end
end


export ChiSquared, Cash
