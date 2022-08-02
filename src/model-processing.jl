export processmodel, freeze!, unfreeze!

struct ProcessedSpectralModel{E,M,P,L}
    expression::E
    models::M
    parameters::P
    modelparams::L
end

function setfreeze!(psm::ProcessedSpectralModel, state::Bool, symbs::Vararg{Symbol})
    foreach(psm.parameters) do (s, p)
        if s in symbs
            p.frozen = state
        end
    end
    psm.parameters
end

function freeze!(psm::ProcessedSpectralModel, symbs::Vararg{Symbol})
    setfreeze!(psm, true, symbs...)
end

function unfreeze!(psm::ProcessedSpectralModel, symbs::Vararg{Symbol})
    setfreeze!(psm, false, symbs...)
end

function modelinfo(m::ProcessedSpectralModel)
    io = IOBuffer()
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == ProcessedSpectralModel
        print(
            io,
            "ProcessedSpectralModel(#Models=$(length(m.models)),#Params=$(length(m.parameters)))",
        )
    else
        write(
            io,
            "ProcessedSpectralModel:\n   $(readable_model_expression(m.expression))\n",
        )

        write(io, "  Models:\n")
        model_infos = map(m.models) do (s, M)
            params = join(m.modelparams[s], ", ")
            Base.typename(M).name, s, params
        end

        model_pad = maximum(length.(String.(first.(model_infos))))
        for (name, s, info) in model_infos
            write(io, "     . $s => $(rpad(name, model_pad))   : $info\n")
        end

        write(io, "  Parameters:\n")
        param_pad = maximum(length.(String.(first.(m.parameters)))) + 1
        for (s, p) in m.parameters
            write(io, "     . $(rpad(s, param_pad)) => $p\n")
        end
    end
    String(take!(io))
end

readable_model_expression(psm::ProcessedSpectralModel) =
    readable_model_expression(psm.expression)
readable_model_expression(symb::Symbol) = symb

function readable_model_expression(expr::Pair)
    recursive_expression_call(expr) do (left, right, op)
        if isnothing(op)
            :($(left)($right))
        else
            if op == :*
                :($left * $right)
            else
                :($left + $right)
            end
        end
    end
end


recursive_expression_call(func, psm::ProcessedSpectralModel) =
    recursive_expression_call(func, psm.expression)

function recursive_expression_call(func, expr::Pair)
    right_expr, (op, left_expr) = expr

    if typeof(right_expr) !== Symbol
        right = recursive_expression_call(func, right_expr)
    else
        right = right_expr
    end
    if typeof(left_expr) !== Symbol
        left = recursive_expression_call(func, left_expr)
    else
        left = left_expr
    end

    func((left, right, op))
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, m::ProcessedSpectralModel)
    print(io, modelinfo(m))
end

# ..... #

function assemble_symbols!(tracker, m::M) where {M<:AbstractSpectralModel}
    if modelkind(M) === Additive
        s = Symbol('a', tracker.additive[])
        push!(tracker.sym_model, s => m)
        tracker.additive[] += 1
        s
    elseif modelkind(M) === Multiplicative
        s = Symbol('m', tracker.multiplicative[])
        push!(tracker.sym_model, s => m)
        tracker.multiplicative[] += 1
        s
    else
        s = Symbol('c', tracker.convolutional[])
        push!(tracker.sym_model, s => m)
        tracker.convolutional[] += 1
        s
    end
end

function assemble_symbols!(tracker, cm::CompositeSpectralModel{M1,M2}) where {M1,M2}
    left = assemble_symbols!(tracker, cm.left)
    right = assemble_symbols!(tracker, cm.right)
    op = if modelkind(M1) === Multiplicative
        :(*)
    elseif modelkind(M1) === Additive
        :(+)
    else
        nothing
    end
    right => (op, left)
end

function addparams!(all_params, m::AbstractSpectralModel)
    map(fieldnames(typeof(m))) do p
        i = 1
        symb = Symbol(p, '_', i)
        while symb in first.(all_params)
            i += 1
            symb = Symbol(p, '_', i)
        end
        push!(all_params, symb => FitParameter(getproperty(m, p)))
        symb
    end
end

function parse_models_with_params!(all_params, models, root_symb::Char)
    map(enumerate(models)) do (i, m)
        params = addparams!(all_params, m)
        instance = Expr(:call, typeof(m), params...)
        Symbol(root_symb, i) => instance
    end
end

function unpack_model_parameters(sym_model)
    all_params = Pair{Symbol,AbstractFitParameter}[]

    model_to_params = map(sym_model) do (s, m)
        model_params = map(fieldnames(typeof(m))) do p
            i = 1
            symb = Symbol(p, '_', i)
            while symb in first.(all_params)
                i += 1
                symb = Symbol(p, '_', i)
            end
            push!(all_params, symb => getproperty(m, p))
            symb
        end
        s => model_params
    end

    Dict(model_to_params), all_params
end

function processmodel(cm::AbstractSpectralModel)
    # assemble evaluation
    tracker = (
        sym_model = Pair{Symbol,AbstractSpectralModel}[],
        additive = Ref(1),
        multiplicative = Ref(1),
        convolutional = Ref(1),
    )

    expr = assemble_symbols!(tracker, cm)
    model_to_params, all_params = unpack_model_parameters(tracker.sym_model)
    models = map(tracker.sym_model) do (s, m)
        s => typeof(m)
    end
    ProcessedSpectralModel(expr, models, all_params, model_to_params)
end

# function build_lsqfit(pm::ProcessedSpectralModel)
#     flux_calls = [:($name = $(model)(energy)) for (name, model) in getmodels(pm)]

#     i = 1
#     params = map(pm.parameters) do (p, f)
#         if f.frozen
#             val = value(f)
#             :($p = $val)
#         else
#             i += 1
#             :($p = params[$(i-1)])
#         end
#     end

#     # params = [:($p = params[$i]) for (i, p) in enumerate(first.(pm.parameters))]
#     func = :((energy, params) -> begin
#         $(params...)
#         $(flux_calls...)
#         res = @. $(pm.expression)
#         @view res[1:end-1]
#     end) |> eval

#     fs = last.(getfreeparams(pm))
#     # model, p0, lower bounds, upper bounds
#     func, value.(fs), lowerbound.(fs), upperbound.(fs)
# end

# export processmodel, ProcessedSpectralModel, getmodels, build_lsqfit, FitParameter, value, freeze!, unfreeze!


# #=
# flux1
# flux2
# flux3

# # m3(m1 * a1 + c1(m2 * a2)) + a3

# a1 = XS_Powerlaw(...)
# #...

# invoke!(flux1, a1)
# invoke!(flux2, m1)
# @. flux1 *= flux2
# invoke!(flux2, a2)
# invoke!(flux3, m2)
# @. flux2 *= flux3
# invoke!(flux2, c1)
# @. flux1 += flux2
# invoke!(flux2, a3)
# @. flux1 += flux3

# return flux1
# =#
