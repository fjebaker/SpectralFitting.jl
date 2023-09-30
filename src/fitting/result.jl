export FittingResult, MultiFittingResult, AbstractFittingResult

abstract type AbstractFittingResult end

struct FittingResult{T,K,C} <: AbstractFittingResult
    χ2::T
    u::K
    config::C
end

function _pretty_print_result(model, u, chi2)
    ppx2 = prettyfloat(chi2)
    ppu = join((prettyfloat(i) for i in u), ", ")
    """
      Model: $(model)
      . u     : [$(ppu)]
      . χ²    : $(ppx2) 
    """
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

struct MultiFittingSlice{C,V,U,T}
    cache::C
    domain::V
    u::U
    χ2::T
end

function Base.getindex(result::MultiFittingResult, i::Int)
    cache = result.config.cache.caches[i]
    u = result.us[i]
    chi2 = result.χ2s[i]
    s, e = _get_range(result.config.cache.domain_mapping, i)
    MultiFittingSlice(cache, result.config.domain[s:e], u, chi2)
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
