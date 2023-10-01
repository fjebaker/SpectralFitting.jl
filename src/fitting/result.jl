export FittingResult,
    MultiFittingResult, AbstractFittingResult, FittingResultSlice, invoke_result

function _pretty_print_result(model, u, chi2)
    ppx2 = prettyfloat(chi2)
    ppu = join((prettyfloat(i) for i in u), ", ")
    """
      Model: $(model)
      . u     : [$(ppu)]
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
    χ2::T
end

measure(stat::AbstractStatistic, slice::FittingResultSlice) = measure(stat, slice, slice.u)

function measure(stat::AbstractStatistic, slice::FittingResultSlice, u)
    measure(stat, slice.objective, invoke_result(slice, u), slice.variance)
end

function invoke_result(slice::FittingResultSlice, u)
    @assert length(u) == length(slice.u)
    _invoke_and_transform!(slice.cache, slice.domain, u)
end

function _pretty_print(slice::FittingResultSlice)
    "FittingResultSlice:\n" * _pretty_print_result(slice.cache.model, slice.u, slice.χ2)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(slice::FittingResultSlice))
    print(io, encapsulate(_pretty_print(slice)))
end

struct FittingResult{T,K,C} <: AbstractFittingResult
    χ2::T
    u::K
    config::C
end

function invoke_result(result::FittingResult, u)
    @assert length(u) == length(result.u)
    _invoke_and_transform!(result.config.cache, result.config.domain, u)
end

function Base.getindex(result::FittingResult, i)
    if i == 1
        FittingResultSlice(
            result.config.cache,
            result.config.domain,
            result.config.objective,
            result.config.variance,
            result.u,
            result.χ2,
        )
    else
        throw(BoundsError())
    end
end

function _pretty_print(res::FittingResult)
    "FittingResult:\n" * _pretty_print_result(res.config.cache.model, res.u, res.χ2)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(res::FittingResult))
    print(io, encapsulate(_pretty_print(res)))
end

struct MultiFittingResult{T,K,C} <: AbstractFittingResult
    χ2s::Vector{T}
    us::K
    config::C
end

function Base.getindex(result::MultiFittingResult, i::Int)
    cache = result.config.cache.caches[i]
    u = result.us[i]
    chi2 = result.χ2s[i]
    d_start, d_end = _get_range(result.config.cache.domain_mapping, i)
    o_start, o_end = _get_range(result.config.cache.objective_mapping, i)
    @views FittingResultSlice(
        cache,
        result.config.domain[d_start:d_end],
        result.config.objective[o_start:o_end],
        result.config.variance[o_start:o_end],
        u,
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
        b = _pretty_print_result(slice.cache.model, slice.u, slice.χ2)
        r = indent(b, 1)
        print(buff, r)
    end
    text = String(take!(buff))
    print(io, encapsulate(text) * "Σχ² = $(total_χ2)")
end
