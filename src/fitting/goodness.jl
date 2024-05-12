
function goodness(
    result::FittingResultSlice,
    u::AbstractVector{T},
    σ::AbstractVector{T};
    N = 1000,
    stat = ChiSquared(),
    distribution = Distributions.Normal,
) where {T}
    x = similar(u)
    measures = zeros(T, N)

    for i in eachindex(measures)
        # sample the next parameters
        for i in eachindex(u)
            m = u[i]
            d = σ[i]
            # TODO: respect the upper and lower bounds of the parameters
            # Distributions.Truncated(distribution(mean, dev), ) 
            x[i] = rand(distribution(m, d))
        end
        measures[i] = measure(stat, result, x)
    end

    perc = 100 * count(<(result.χ2), measures) / N
    @info "% with measure < result = $(perc)"

    measures
end 

function goodness(result::FittingResult, σu = result.σu; kwargs...)
    @assert !isnothing(result.σu) "σ cannot be nothing, else algorithm has no parameter intervals to sample from."
    goodness(result[1], result.u, σu; kwargs...)
end

export goodness