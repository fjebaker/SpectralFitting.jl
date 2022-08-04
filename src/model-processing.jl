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

get_params(psm::ProcessedSpectralModel, symb...) = map(i -> get_params(psm, i), symb)

function set_freeze!(psm::ProcessedSpectralModel, state::Bool, symbs::Vararg{Symbol})
    action! = state ? freeze! : unfreeze!
    foreach(psm.parameters) do (s, p)
        if s in symbs
            action!(p)
        end
    end
    psm.parameters
end

freeze!(psm::ProcessedSpectralModel, symbs::Vararg{Symbol}) =
    setfreeze!(psm, true, symbs...)
unfreeze!(psm::ProcessedSpectralModel, symbs::Vararg{Symbol}) =
    setfreeze!(psm, false, symbs...)
Base.show(io::IO, ::MIME{Symbol("text/plain")}, m::ProcessedSpectralModel) =
    print_info(io, m)
Base.show(io::IO, m::ProcessedSpectralModel) = print(io, "ProcessedSpectralModel")

function print_info(io::IO, m::ProcessedSpectralModel)
    header = "ProcessedSpectralModel:\n   $(readable_model_expression(m.expression))\n"
    print(io, header)

    print(io, "  Models:\n")
    model_infos = map(m.models) do (s, M)
        params = join(m.modelparams[s], ", ")
        String(model_type_name(M)), s, params
    end

    model_pad = maximum(length.(first.(model_infos)))
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

#=
    Parsing Utilities
=#
recursive_expression_call(func, psm::ProcessedSpectralModel) =
    recursive_expression_call(func, psm.expression)

function recursive_expression_call(func, expr::Pair)
    right_expr, (op, left_expr) = expr
    right =
        typeof(right_expr) !== Symbol ? recursive_expression_call(func, right_expr) :
        right_expr
    left =
        typeof(left_expr) !== Symbol ? recursive_expression_call(func, left_expr) :
        left_expr
    func((left, right, op))
end

function add_model_symbol!(c::Char, ref)
    i = (ref[] += 1)
    Symbol(c, i)
end
process_model_symbols!(tracker, ::M) where {M} =
    process_model_symbols!(tracker, modelkind(M))
process_model_symbols!(tracker, ::Additive) = add_model_symbol!('a', tracker.additive)
process_model_symbols!(tracker, ::Multiplicative) =
    add_model_symbol!('m', tracker.multiplicative)
process_model_symbols!(tracker, ::Convolutional) =
    add_model_symbol!('c', tracker.convolutional)

function assemble_symbols!(tracker, m::M) where {M<:AbstractSpectralModel}
    s = process_model_symbols!(tracker, m)
    push!(tracker.model_symbols, s => m)
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

function make_unique_param_symbol(p, all_symbols)
    i = 1
    symb = Symbol(p, '_', i)
    while symb in all_symbols
        i += 1
        symb = Symbol(p, '_', i)
    end
    symb
end

function unpack_model_parameters(model_symbols, T::Type)
    all_parameters = Pair{Symbol,FitParam{T}}[]
    model_symbol_to_parameters = map(model_symbols) do (model_symbol, m)
        model_params = map(parameter_symbol_pairs(m)) do (param_symbol, param)
            symb = make_unique_param_symbol(param_symbol, first.(all_parameters))
            push!(all_parameters, symb => deepcopy(param))
            symb
        end
        model_symbol => model_params
    end
    Dict(model_symbol_to_parameters), all_parameters
end

function processmodel(cm::AbstractSpectralModel)
    # assemble evaluation
    tracker = (
        model_symbols = Pair{Symbol,AbstractSpectralModel}[],
        additive = Ref(0),
        multiplicative = Ref(0),
        convolutional = Ref(0),
    )
    expr = assemble_symbols!(tracker, cm)
    model_to_params, all_params =
        unpack_model_parameters(tracker.model_symbols, numbertype(cm))
    models = map(tracker.model_symbols) do (s, m)
        s => typeof(m)
    end
    ProcessedSpectralModel(expr, models, all_params, model_to_params)
end
