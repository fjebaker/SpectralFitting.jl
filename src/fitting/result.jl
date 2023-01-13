
struct FittingResult{T,M,E,F}
    u::Vector{T}
    χ2::T
    model::M
    x::E
    folded_invoke::F
end

function Base.show(io::IO, ::MIME"text/plain", res::FittingResult)
    println(
        io,
        """FittingResult:
          Model: $(res.model)
            - Parameters  : $(res.u)
            - χ2          : $(res.χ2) 
        """,
    )
end

function bundle_result(u, model, f, x, y, variance)
    chi2 = χ2_from_ŷyvar(f(x, u), y, variance)
    FittingResult(u, chi2, model, x, f)
end

struct MultiFittingResult{F}
    results::F
    MultiFittingResult(result::Tuple) = new{typeof(result)}(result)
end

function Base.show(io::IO, ::MIME"text/plain", res::MultiFittingResult)
    println(io, "┌ MultiFittingResult:")
    buff = IOBuffer()
    for (i, result) in enumerate(res.results)
        show(buff, MIME"text/plain"(), result)
        s = String(take!(buff))
        s = s[1:findlast(==('\n'), s)-1]
        text = if i != lastindex(res.results)
            replace(s, "\n" => "\n│ ")
        else
            _s = replace(s, "\n" => "\n│ ")
            _s[1:findlast(==('│'), _s)-1] * '└'
        end
        println(io, "│ " * text)
    end
end

function unpack_multimodel(parameters, m::MultiModel, X, Y, V, state)
    n_params = state.i_params
    n_energy = state.i_x
    n_output = state.i_out
    results = map((1:state.n_models...,)) do i
        model = m.m[i]
        f = state.funcs[i]
        start_p = i == 1 ? 1 : n_params[i-1] + 1
        end_p = n_params[i]
        # don't view here as we want a copy for the output
        u = parameters[start_p:end_p]

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
