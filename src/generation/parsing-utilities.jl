const Lens = Union{Symbol,Expr}

struct ModelInfo
    symbols::Vector{Symbol}
    free::Vector{Symbol}
    frozen::Vector{Symbol}
    generated_symbols::Vector{Symbol}
    lens::Lens
    type::Type
end

mutable struct GenerationAggregate{NumType}
    statements::Vector{Expr}
    infos::Vector{ModelInfo}
    closure_params::Vector{Symbol}
    models::Vector{Type}
    flux_count::Int
    maximum_flux_count::Int
end
GenerationAggregate(T) =
    GenerationAggregate{T}(Expr[], ModelInfo[], Symbol[], Type[], 0, 0)
GenerationAggregate(T) =
    GenerationAggregate{T}(Expr[], ModelInfo[], Symbol[], Type[], 0, 0)

push_closure_param!(g::GenerationAggregate, s::Symbol) = push!(g.closure_params, s)
push_info!(g::GenerationAggregate, s::ModelInfo) = push!(g.infos, s)
push_statement!(g::GenerationAggregate, s::Expr) = push!(g.statements, s)
push_model!(g::GenerationAggregate, t::Type) = push!(g.models, t)

function set_flux!(g::GenerationAggregate, f)
    g.flux_count = f
    if g.flux_count > g.maximum_flux_count
        g.maximum_flux_count = g.flux_count
    end
    g.flux_count
end

inc_flux!(g::GenerationAggregate) = set_flux!(g, g.flux_count + 1)
dec_flux!(g::GenerationAggregate) = set_flux!(g, g.flux_count - 1)

get_flux_symbol(i::Int) = Symbol(:flux, i)

function getinfo(model::Type{<:AbstractSpectralModel}; lens::Lens = :(model))
    symbs = [all_parameter_symbols(model)...]
    free = [free_parameter_symbols(model)...]
    frozen = [frozen_parameter_symbols(model)...]
    ModelInfo(symbs, free, frozen, [Base.gensym(s) for s in symbs], lens, model)
end

function getinfo(model::Type{<:CompositeModel}; lens::Lens = :(model))
    infos = ModelInfo[]
    _addinfo!(infos, model, lens)
    infos
end

function _addinfo!(infos, model::Type{<:AbstractSpectralModel}, lens::Lens)
    push!(infos, getinfo(model; lens = lens))
end
function _addinfo!(infos, model::Type{<:CompositeModel}, lens::Lens)
    right_run = :(getproperty($lens, :right))
    _addinfo!(infos, model.parameters[2], right_run)
    left_run = :(getproperty($lens, :left))
    _addinfo!(infos, model.parameters[1], left_run)
end

function _addinfosymbol!(infos, counters, model::Type{<:AbstractSpectralModel}, lens::Lens)
    symb = _unique_model_symbol(modelkind(model), counters)
    push!(infos, symb => getinfo(model; lens = lens))
    symb
end
function _addinfosymbol!(infos, counters, model::Type{<:CompositeModel}, lens::Lens)
    right_run = :(getproperty($lens, :right))
    r_symb = _addinfosymbol!(infos, counters, model.parameters[2], right_run)
    left_run = :(getproperty($lens, :left))
    l_symb = _addinfosymbol!(infos, counters, model.parameters[1], left_run)
    op = operation_symbol(model.parameters[3])
    if isnothing(op)
        :($(l_symb)($r_symb))
    else
        Expr(:call, op, l_symb, r_symb)
    end
end

function _addinfoinvoke!(ga::GenerationAggregate{NumType}, model::Type{<:AbstractSpectralModel}, lens::Lens) where {NumType}
    # don't increment flux for convolutional models
    if !(modelkind(model) === Convolutional)
        inc_flux!(ga)
    end
    flux = get_flux_symbol(ga.flux_count)
    # get and append model info
    info = getinfo(model; lens = lens)
    push!(ga.infos, info)

    # aggregate closure parameters
    closure_params = if has_closure_params(model)
        map(closure_parameter_symbols(model)) do p
            new_closure_param!(ga, p)
        end
    else
        ()
    end

    # get the parameter type
    T = NumType <: SpectralFitting.FitParam ? NumType.parameters[1] : NumType
    model_constructor =
        :($(model.name.wrapper){$(model.parameters[1:end-2]...),$(T),$(model.parameters[end])})

    # assemble the invocation statement
    s = :(invokemodel!(
        $flux,
        energy,
        $(model_constructor)($(closure_params...), $(info.generated_symbols...)),
    ))
    push_model!(ga, model)
    push_statement!(ga, s)
end

function _addinfoinvoke!(ga::GenerationAggregate, model::Type{<:CompositeModel}, lens::Lens)
    right_run = :(getproperty($lens, :right))
    _addinfoinvoke!(ga, model.parameters[2], right_run)
    left_run = :(getproperty($lens, :left))
    _addinfoinvoke!(ga, model.parameters[1], left_run)

    # resolve the flux pair using the reduction operation of the models
    op = operation_symbol(model.parameters[3])
    if !isnothing(op)
        add_flux_resolution!(ga, op)
    end
end

function add_flux_resolution!(ga::GenerationAggregate, op::Symbol)
    fr = get_flux_symbol(ga.flux_count)
    dec_flux!(ga)
    fl = get_flux_symbol(ga.flux_count)
    expr = Expr(:call, op, fl, fr)
    push_statement!(ga, :(@.($fl = $expr)))
end

function _make_unique_readable_symbol(p, symbol_bank; delim = '_')
    i = 1
    symb = Symbol(p, delim, i)
    while symb in symbol_bank
        i += 1
        symb = Symbol(p, delim, i)
    end
    symb
end

function all_parameter_symbols(M::Type{<:AbstractSpectralModel})
    fieldnames(M)
end

function free_parameter_symbols(M::Type{<:AbstractSpectralModel})
    first(M.parameters[end].parameters)
end

function frozen_parameter_symbols(M::Type{<:AbstractSpectralModel})
    free_params = free_parameter_symbols(M)
    frozen = filter(i -> i ∉ free_params, all_parameter_symbols(M))
    tuple(frozen...)
end

all_parameter_symbols(::Type{<:CompositeModel}) = error("Unreachable.")
free_parameter_symbols(::Type{<:CompositeModel}) = error("Unreachable.")
frozen_parameter_symbols(::Type{<:CompositeModel}) = error("Unreachable.")

function closure_parameter_symbols(::Type{<:AbstractSpectralModel})
    error("This specialisation should never need to be invoked.")
end

model_base_name(M::Type{<:AbstractSpectralModel}) = Base.typename(M).name

function _vector_to_named_tuple(params, names)
    statements = [:(params[$(i)]) for i = 1:length(names)]
    :(NamedTuple{$(names)}(($(statements...),)))
end

function free_parameters_to_named_tuple(params, model)
    names = free_parameter_symbols(model)
    _vector_to_named_tuple(params, names)
end

function frozen_parameters_to_named_tuple(params, model)
    names = frozen_parameter_symbols(model)
    _vector_to_named_tuple(params, names)
end

function all_parameters_to_named_tuple(params, model::Type{<:AbstractSpectralModel})
    names = all_parameter_symbols(model)
    _vector_to_named_tuple(params, names)
end

function all_parameters_to_named_tuple(model::Type{<:AbstractSpectralModel})
    names = all_parameter_symbols(model)
    statements = [:(getproperty(model, $(Meta.quot(s)))) for s in names]
    :(NamedTuple{$(names)}(($(statements...),)))
end

function _unique_parameter_symbols(infos::Vector{ModelInfo})
    all_parameters = reduce(vcat, map(i -> i.symbols, infos))
    param_names = Symbol[]
    for param in all_parameters
        s = _make_unique_readable_symbol(param, param_names)
        push!(param_names, s)
    end
    (param_names...,)
end

function all_parameters_to_named_tuple(model::Type{<:CompositeModel})
    infos = getinfo(model)
    lenses = reduce(vcat, map(i -> _parameter_lens(i, i.symbols), infos))
    names = _unique_parameter_symbols(infos)
    :(NamedTuple{$(names)}(($(lenses...),)))
end

_unique_model_symbol(::SpectralFitting.Additive, counters) = Symbol('a', counters.a[] += 1)
_unique_model_symbol(::SpectralFitting.Multiplicative, counters) =
    Symbol('m', counters.m[] += 1)
_unique_model_symbol(::SpectralFitting.Convolutional, counters) =
    Symbol('c', counters.c[] += 1)

function _destructure_for_printing(model::Type{<:CompositeModel}; lens = :(model))
    counters = (; a = Ref(0), m = Ref(0), c = Ref(0))
    infos = Pair{Symbol,ModelInfo}[]
    expr = _addinfosymbol!(infos, counters, model, lens)
    # reorder data structure as NamedTuple of pairs of (model, parameters)
    param_names = _unique_parameter_symbols(last.(infos))
    offset = 1
    descriptions = map(infos) do (name, info)
        # get parameter symbols
        params = all_parameter_symbols(info.type)
        free_params = free_parameter_symbols(info.type)

        # get the symbols generated for the composite model
        n_params = length(params) - 1
        names = param_names[offset:offset+n_params]
        offset += n_params + 1

        descr = :(
            $(info.lens),
            ($(map(i -> "$i", names)...),),
            ($(map(i -> i in free_params, params)...),),
        )
        descr
    end

    expr_string = "$expr"

    quote
        k = NamedTuple{($(Meta.quot.(first.(infos))...),)}(($(descriptions...),))
        $(expr_string), k
    end
end
