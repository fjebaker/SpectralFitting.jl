export processmodel, freeze!, unfreeze!, setparams!

struct ProcessedSpectralModel{E,M,P,L}
    expression::E
    models::M
    parameters::P
    modelparams::L
end

function get_params(psm::ProcessedSpectralModel, symb)
    i = findfirst(i -> first(i) == symb, psm.parameters)
    if !isnothing(i)
        psm.parameters[i]
    else
        error("No parameter with symbol: $symb.")
    end
end

get_params(psm::ProcessedSpectralModel, symb...) =
    map(i -> get_params(psm, i), symb)

function set_freeze!(psm::ProcessedSpectralModel, state::Bool, symbs::Vararg{Symbol})
    action! = state ? freeze! : unfreeze!
    foreach(psm.parameters) do (s, p)
        if s in symbs
            action!(p)
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

function print_info(io::IO, m::ProcessedSpectralModel)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == ProcessedSpectralModel
        print(
            io,
            "ProcessedSpectralModel(#Models=$(length(m.models)),#Params=$(length(m.parameters)))",
        )
    else
        print(
            io,
            "ProcessedSpectralModel:\n   $(readable_model_expression(m.expression))\n",
        )

        print(io, "  Models:\n")
        model_infos = map(m.models) do (s, M)
            params = join(m.modelparams[s], ", ")
            Base.typename(M).name, s, params
        end

        model_pad = maximum(length.(String.(first.(model_infos))))
        for (name, s, info) in model_infos
            print(io, "     . $s => $(rpad(name, model_pad))   : $info\n")
        end

        print(io, "  Parameters:\n")
        param_symbs = String.(first.(m.parameters))
        param_instance = last.(m.parameters)
        param_pad = maximum(length.(param_symbs)) + 1
        param_infos = get_info_tuple.(param_instance)
        frozen_states = is_frozen.(param_instance)
        pads = map(1:2) do j
            maximum(i -> length(i[j]), param_infos) + 1
        end
        for (s, p, f) in zip(param_symbs, param_infos, frozen_states)
            print(io, "     . $(rpad(s, param_pad)) => ")
            print(io, rpad(p[1], pads[1]), " ∈ ", rpad(p[2], pads[2]))
            state = if f
                Crayons.Box.CYAN_FG("Frozen")
            else
                Crayons.Box.GREEN_FG("Free")
            end
            print(io, state)
            print(io, "\n")
        end
    end
    io
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
    print_info(io, m)
end

# ..... #

process_model_symbols!(tracker, ::M) where {M} = process_model_symbols!(tracker, modelkind(M))
function process_model_symbols!(tracker, ::Additive)
    s = Symbol('a', tracker.additive[])
    tracker.additive[] += 1
    s
end
function process_model_symbols!(tracker, ::Multiplicative)
    s = Symbol('m', tracker.multiplicative[])
    tracker.multiplicative[] += 1
    s
end
function process_model_symbols!(tracker, ::Convolutional)
    s = Symbol('c', tracker.convolutional[])
    tracker.convolutional[] += 1
    s
end

function assemble_symbols!(tracker, m::M) where {M<:AbstractSpectralModel}
    s = process_model_symbols!(tracker, m)
    push!(tracker.sym_model, s => m)
    s
end

get_operation(::Multiplicative) = :(*)
get_operation(::Additive) = :(+)
get_operation(::Convolutional) = nothing

function assemble_symbols!(tracker, cm::CompositeSpectralModel{M1,M2}) where {M1,M2}
    left = assemble_symbols!(tracker, cm.left)
    right = assemble_symbols!(tracker, cm.right)
    op = get_operation(modelkind(M1))
    right => (op, left)
end

function unpack_model_parameters(sym_model, T::Type)
    all_params = Pair{Symbol,FitParam{T}}[]

    model_to_params = map(sym_model) do (s, m)
        model_params = map(fieldnames(typeof(m))) do p
            i = 1
            symb = Symbol(p, '_', i)
            while symb in first.(all_params)
                i += 1
                symb = Symbol(p, '_', i)
            end
            push!(all_params, symb => deepcopy(getproperty(m, p)))
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
    model_to_params, all_params = unpack_model_parameters(tracker.sym_model, numbertype(cm))
    models = map(tracker.sym_model) do (s, m)
        s => typeof(m)
    end
    ProcessedSpectralModel(expr, models, all_params, model_to_params)
end
