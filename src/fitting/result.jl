export FittingResult,
    MultiFittingResult,
    AbstractFittingResult,
    FittingResultSlice,
    invoke_result,
    update_model!

function _pretty_print_result(model, u, σ, chi2)
    ppx2 = prettyfloat(chi2)
    ppu = join((prettyfloat(i) for i in u), ", ")
    ppσ = isnothing(σ) ? nothing : join((prettyfloat(i) for i in σ), ", ")
    """
      Model: $(model)
      . u     : [$(ppu)]
      . σᵤ    : [$(ppσ)]
      . χ²    : $(ppx2) 
    """
end

abstract type AbstractFittingResult end

invoke_result(res::AbstractFittingResult) = invoke_result(res, res.u)

struct FittingResultSlice{C,V,U,T} <: AbstractFittingResult
    cache::C
    domain::V
    objective::V
    variance::V
    u::U
    σu::Union{Nothing,U}
    χ2::T
end

estimated_error(r::FittingResultSlice) = r.σu
estimated_params(r::FittingResultSlice) = r.u
measure(stat::AbstractStatistic, slice::FittingResultSlice) = measure(stat, slice, slice.u)

function measure(stat::AbstractStatistic, slice::FittingResultSlice, u)
    measure(stat, slice.objective, invoke_result(slice, u), slice.variance)
end

function invoke_result(slice::FittingResultSlice, u)
    @assert length(u) == length(slice.u)
    _invoke_and_transform!(slice.cache, slice.domain, u)
end

function _pretty_print(slice::FittingResultSlice)
    "FittingResultSlice:\n" *
    _pretty_print_result(slice.cache.model, slice.u, slice.σu, slice.χ2)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(slice::FittingResultSlice))
    print(io, encapsulate(_pretty_print(slice)))
end

struct FittingResult{T,U,C} <: AbstractFittingResult
    χ2::T
    u::U
    σu::Union{Nothing,U}
    config::C
end

estimated_error(r::FittingResult) = r.σu
estimated_params(r::FittingResult) = r.u

measure(stat::AbstractStatistic, slice::FittingResult, args...) =
    measure(stat, slice[1], args...)

function invoke_result(result::FittingResult, u)
    @assert length(u) == length(result.u)
    _invoke_and_transform!(result.config.cache, result.config.model_domain, u)
end

function Base.getindex(result::FittingResult, i)
    if i == 1
        FittingResultSlice(
            result.config.cache,
            result.config.model_domain,
            result.config.objective,
            result.config.variance,
            result.u,
            result.σu,
            result.χ2,
        )
    else
        throw(BoundsError())
    end
end

function _pretty_print(res::FittingResult)
    "FittingResult:\n" * _pretty_print_result(res.config.cache.model, res.u, res.σu, res.χ2)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(res::FittingResult))
    print(io, encapsulate(_pretty_print(res)))
end

struct MultiFittingResult{T,U,C} <: AbstractFittingResult
    χ2s::Vector{T}
    us::U
    σus::Union{Nothing,U}
    config::C
end

estimated_error(r::MultiFittingResult) = r.σus
estimated_params(r::MultiFittingResult) = r.us

function Base.getindex(result::MultiFittingResult, i::Int)
    cache = result.config.cache.caches[i]
    u = result.us[i]
    σu = isnothing(result.σus) ? nothing : result.σus[i]
    chi2 = result.χ2s[i]
    d_start, d_end = _get_range(result.config.cache.domain_mapping, i)
    o_start, o_end = _get_range(result.config.cache.objective_mapping, i)
    @views FittingResultSlice(
        cache,
        result.config.model_domain[d_start:d_end],
        result.config.objective[o_start:o_end],
        result.config.variance[o_start:o_end],
        u,
        σu,
        chi2,
    )
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(res::MultiFittingResult))
    total_χ2 = prettyfloat(sum(res.χ2s))

    buff = IOBuffer()
    println(buff, "MultiFittingResult:")
    print(buff, " ")
    for i = 1:length(res.us)
        slice = res[i]
        b = _pretty_print_result(slice.cache.model, slice.u, slice.σu, slice.χ2)
        r = indent(b, 1)
        print(buff, r)
    end
    text = String(take!(buff))
    print(io, encapsulate(text) * "Σχ² = $(total_χ2)")
end

function update_model!(
    model::AbstractSpectralModel,
    result::Union{<:FittingResult,<:FittingResultSlice},
)
    free_params = filter(isfree, parameter_tuple(model))
    for (i, f) in enumerate(free_params)
        set_value!(f, result.u[i])
    end
    model
end

function update_model!(
    model::FittableMultiModel,
    result::Union{<:FittingResult,<:FittingResultSlice},
)
    @assert length(model.m) == 1
    update_model!(model.m[1], result)
end

function update_model!(multimodel::FittableMultiModel, result::MultiFittingResult)
    error("TODO")
end
