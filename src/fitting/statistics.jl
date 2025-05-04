
function measure(::ChiSquared, y, ŷ, σ²)
    sum(@.((y - ŷ)^2 / σ²))
end

function measure(::Cash, S, ŷ, _) # ignores variance
    # avoid errors in log
    any(<(0), ŷ) && return Inf
    any(<(0), S) && return Inf
    2 * sum(@.(ŷ - S + S * (log(S) - log(ŷ))))
end

function measure(s::AbstractStatistic, slice::FitResultSlice, u = slice.u)
    ŷ = calculate_objective!(slice, u)
    measure(s, get_objective(slice), ŷ, get_objective_variance(slice))
end

function measure(stat::AbstractStatistic, result::FitResult, args...; kwargs...)
    measure(stat, result[1], args...; kwargs...)
end

export ChiSquared, Cash, measure
