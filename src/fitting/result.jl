export FittingResult, MultiFittingResult, AbstractFittingResult

abstract type AbstractFittingResult end

struct FittingResult{T,M,E,F} <: AbstractFittingResult
    u::Vector{T}
    χ2::T
    model::M
    x::E
    folded_invoke::F
end

function _pretty_print(res::FittingResult)
    chi2 = prettyfloat(res.χ2)
    us = join((prettyfloat(i) for i in res.u), ", ")
    """FittingResult:
        Model: $(res.model)
        . u     : [$(us)]
        . χ²    : $(chi2) 
    """
end

function Base.show(io::IO, ::MIME"text/plain", res::FittingResult)
    print(io, encapsulate(_pretty_print(res)))
end

function bundle_result(u, model, f, x, y, variance)
    chi2 = χ2_from_ŷyvar(f(x, u), y, variance)
    FittingResult(u, chi2, model, x, f)
end

struct MultiFittingResult{F} <: AbstractFittingResult
    results::F
    MultiFittingResult(result::Tuple) = new{typeof(result)}(result)
end

function Base.show(io::IO, ::MIME"text/plain", res::MultiFittingResult)
    total_χ2 = prettyfloat(sum(i -> i.χ2, res.results))

    buff = IOBuffer()
    println(buff, "MultiFittingResult:")
    print(buff, " ")
    for result in res.results
        b = _pretty_print(result) * "\n"
        r = indent(b, 1)
        # drop last new line

        print(buff, r)
    end
    text = String(take!(buff))
    print(io, encapsulate(text) * "Σχ² = $(total_χ2)")
end

function unpack_multimodel(parameters, m::MultiModel, X, Y, V, state)
    parameter_indices = state.parameter_indices
    n_energy = state.i_x
    n_output = state.i_out
    results = map((1:state.n_models...,)) do i
        model = m.m[i]
        f = state.funcs[i]
        # don't view here as we want a copy for the output
        u = parameters[parameter_indices[i]]

        start_x = i == 1 ? 1 : n_energy[i-1] + 1
        end_x = n_energy[i]
        x = X[start_x:end_x]

        start_y = i == 1 ? 1 : n_output[i-1] + 1
        end_y = n_output[i]
        y = Y[start_y:end_y]
        variance = V[start_y:end_y]
        bundle_result(u, model, f, x, y, variance)
    end
    MultiFittingResult(results)
end
