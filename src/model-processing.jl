export processmodel, freeze!, unfreeze!, setparams!

struct ProcessedSpectralModel{E,M,P,L,I}
    expression::E
    model_type_index::M
    parameters::P
    modelparams::L
    instances::I
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

    print(io, "  Models:")
    if !isnothing(m.instances)
        N_non_triv = length(m.instances)
        print(io, "  ($N_non_triv non trivial)")
    end
    print(io, "\n")
    model_infos = map(m.model_type_index) do (s, M)
        params = join(m.modelparams[s], ", ")
        String(model_type_name(M)), s, params
    end

    model_pad = maximum(length.(first.(model_infos)))
    info_pad = maximum(i -> length(i[3]), model_infos) + 2
    for (name, s, info) in model_infos
        print(io, "     . $s => $(rpad(name, model_pad))   : $(rpad(info,info_pad))")
        if (!isnothing(m.instances)) && s in first.(m.instances)
            print(io, Crayons.Box.YELLOW_FG("NON TRIVIAL"))
        end
        print(io, "\n")
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
    left_expr, (op, right_expr) = expr
    # very important is order here:
    # left must be unpacked first to parse multiplicatives before additives
    left =
        typeof(left_expr) !== Symbol ? recursive_expression_call(func, left_expr) :
        left_expr
    right =
        typeof(right_expr) !== Symbol ? recursive_expression_call(func, right_expr) :
        right_expr
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

function assemble_expression!(tracker, m::M) where {M<:AbstractSpectralModel}
    s = process_model_symbols!(tracker, m)
    push!(tracker.symbol_to_model, s => m)
    s
end

get_operation(::Multiplicative) = :(*)
get_operation(::Additive) = :(+)
get_operation(::Convolutional) = nothing

function assemble_expression!(tracker, cm::CompositeSpectralModel{M1,M2}) where {M1,M2}
    left = assemble_expression!(tracker, cm.left)
    right = assemble_expression!(tracker, cm.right)
    op = get_operation(modelkind(M1))
    left => (op, right)
end

function make_unique_symbol(p, symbol_bank; delim = '_')
    i = 1
    symb = Symbol(p, delim, i)
    while symb in symbol_bank
        i += 1
        symb = Symbol(p, delim, i)
    end
    symb
end

function unpack_model_parameters(symbol_to_model, T::Type)
    all_parameters = Pair{Symbol,FitParam{T}}[]
    model_symbol_to_parameters = map(symbol_to_model) do (model_symbol, m)
        model_params = map(parameter_symbol_pairs(m)) do (param_symbol, param)
            symb = make_unique_symbol(param_symbol, first.(all_parameters))
            push!(all_parameters, symb => deepcopy(param))
            symb
        end
        model_symbol => model_params
    end
    Dict(model_symbol_to_parameters), all_parameters
end

function processmodel(cm::AbstractSpectralModel; deepcopy_nontrivial = false)
    # assemble evaluation
    tracker = (
        symbol_to_model = Pair{Symbol,AbstractSpectralModel}[],
        additive = Ref(0),
        multiplicative = Ref(0),
        convolutional = Ref(0),
    )
    expr = assemble_expression!(tracker, cm)
    model_to_params, all_params =
        unpack_model_parameters(tracker.symbol_to_model, numbertype(cm))
    symbol_to_model_types = map(tracker.symbol_to_model) do (s, m)
        s => typeof(m)
    end
    non_trivial_instances =
        filter(!is_trivially_constructed ∘ last, tracker.symbol_to_model)
    if deepcopy_nontrivial
        non_trivial_instances = deepcopy.(non_trivial_instances)
    end
    ProcessedSpectralModel(
        expr,
        symbol_to_model_types,
        all_params,
        model_to_params,
        non_trivial_instances,
    )
end
