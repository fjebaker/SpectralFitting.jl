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

function Base.getindex(result::MultiFittingResult, i::Int)
    model = result.config.cache.caches[i].model
    u = result.us[i]
    chi2 = result.χ2s[i]
    return (; model = model, u = u, χ2 = chi2)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(res::MultiFittingResult))
    total_χ2 = prettyfloat(sum(res.χ2s))

    buff = IOBuffer()
    println(buff, "MultiFittingResult:")
    print(buff, " ")
    for i = 1:length(res.us)
        model, u, chi2 = res[i]
        b = _pretty_print_result(model, u, chi2)
        r = indent(b, 1)
        print(buff, r)
    end
    text = String(take!(buff))
    print(io, encapsulate(text) * "Σχ² = $(total_χ2)")
end
