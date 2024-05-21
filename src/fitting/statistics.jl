function measure(stat::AbstractStatistic, result::FittingResult, args...; kwargs...)
    measure(stat, result[1], args...; kwargs...)
end

struct ChiSquared <: AbstractStatistic end
function measure(::ChiSquared, y, ŷ, σ²)
    sum(@.((y - ŷ)^2 / σ²))
end
function measure(s::ChiSquared, ::FittingConfig, y, ŷ, σ²)
    measure(s, y, ŷ, σ²)
end
function measure(s::ChiSquared, slice::FittingResultSlice, u = slice.u)
    ŷ = _invoke_and_transform!(get_cache(slice), slice.domain, u)
    measure(s, slice.objective, ŷ, slice.variance)
end

function _f_wrap_objective(stat::ChiSquared, config::FittingConfig)
    f = _f_objective(config)
    function _objective(parameters, domain)
        ŷ = f(domain, parameters)
        measure(stat, config.objective, ŷ, config.variance)
    end
end

struct Cash <: AbstractStatistic end
function measure(::Cash, S, ŷ)
    2 * sum(@.(ŷ - S + S * (log(S) - log(ŷ))))
end
function measure(s::Cash, config::FittingConfig, y, ŷ, σ²)
    exp_time = config.prob.data[1].data.spectrum.exposure_time
    prefactor = exp_time .* bin_widths(config.prob.data[1].data)
    measure(s, y .* prefactor, ŷ, σ²)
end

function _f_wrap_objective(stat::Cash, config::FittingConfig)
    f = _f_objective(config)
    exp_time = config.prob.data[1].spectrum.exposure_time
    prefactor = exp_time .* bin_widths(config.prob.data[1])
    function _objective(parameters, domain)
        ŷ = f(domain, parameters)
        measure(stat, config.objective .* prefactor, ŷ .* prefactor)
    end
end

function measure(s::Cash, slice::FittingResultSlice, u = slice.u)
    data = get_dataset(slice)
    exp_time = data.data.spectrum.exposure_time
    prefactor = exp_time .* bin_widths(data.data)
    ŷ = _invoke_and_transform!(get_cache(slice), slice.domain, u)
    measure(s, slice.objective .* prefactor, y .* prefactor̂)
end

export ChiSquared, Cash, measure
