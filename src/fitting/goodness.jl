function goodness(
    result::FitResult,
    u::AbstractVector{T},
    σ::AbstractVector{T};
    N = 1000,
    stat = ChiSquared(),
    distribution = Distributions.Normal,
    refit = true,
    seed = abs(randint()),
    kwargs...,
) where {T}
    x = similar(u)
    measures = zeros(T, N)
    config = deepcopy(result.config)

    rng = Random.default_rng()
    Random.seed!(rng, seed)

    for i in eachindex(measures)
        # sample the next parameters
        for i in eachindex(u)
            m = u[i]
            d = σ[i]
            # TODO: respect the upper and lower bounds of the parameters
            distr = Distributions.Truncated(
                distribution(m, d),
                get_lowerlimit(config.parameters[i]),
                get_upperlimit(config.parameters[i]),
            )
            x[i] = rand(distr)
        end

        simulate!(config, x; rng = rng, kwargs...)

        if refit
            new_result = fit(config, LevenbergMarquadt())
            measures[i] = measure(stat, new_result)
        else
            measures[i] =
                measure(stat, config.objective, invoke_result(result, x), config.variance)
        end
    end

    perc = 100 * count(<(result.χ2), measures) / N
    @info "% with measure < result = $(perc)"

    measures
end

function goodness(result::FitResult, σu = estimated_error(result); kwargs...)
    @assert !isnothing(σu) "σ cannot be nothing, else algorithm has no parameter intervals to sample from."
    goodness(result, estimated_params(result), σu; kwargs...)
end

export goodness
