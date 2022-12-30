struct ModelInfo
    symbols::Vector{Symbol}
    free::Vector{Symbol}
    frozen::Vector{Symbol}
    generated_symbols::Vector{Symbol}
end

function getinfo(model::Type{<:AbstractSpectralModel})
    symbs = [all_parameter_symbols(model)...]
    free = [free_parameter_symbols(model)...]
    frozen = [frozen_parameter_symbols(model)...]
    ModelInfo(symbs, free, frozen, [Base.gensym(s) for s in symbs])
end

function getinfo(model::Type{<:CompositeModel})
    infos = ModelInfo[]
    addinfo!(infos, model)
    infos
end

function addinfo!(infos, model::Type{<:AbstractSpectralModel})
    push!(infos, getinfo(model))
end
function addinfo!(infos, model::Type{<:CompositeModel})
    FunctionGeneration.recursive_model_parse(model) do (left, right, _)
        if (right !== Nothing)
            addinfo!(infos, right)
        end
        if (left !== Nothing)
            addinfo!(infos, left)
        end
        Nothing
    end
end

unpack_model(::Type{<:CompositeModel{M1,M2,O}}) where {M1,M2,O} = (M1, M2, O)
unpack_model(m::CompositeModel{M1,M2,O}) where {M1,M2,O} =
    m.left, m.right, operation_symbol(O)

@inline function __recursive_model_parse(callback, model)
    left_model, right_model, operator = unpack_model(model)
    left = recursive_model_parse(callback, left_model)
    right = recursive_model_parse(callback, right_model)
    callback((left, right, operator))
end

recursive_model_parse(_, model::Type{<:AbstractSpectralModel}) = model
function recursive_model_parse(callback, model::Type{<:CompositeModel})
    __recursive_model_parse(callback, model)
    # Nothing
end

recursive_model_parse(_, model::AbstractSpectralModel) = model
function recursive_model_parse(callback, model::CompositeModel)
    __recursive_model_parse(callback, model)
    # nothing
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

function index_models!(index, running, ::Type{<:CompositeModel{M1,M2}}) where {M1,M2}
    left_run = :(getproperty($running, :left))
    right_run = :(getproperty($running, :right))
    # order is very imporant!!
    # have to recursively go down left branch then right
    if M1 <: CompositeModel
        index_models!(index, left_run, M1)
    end
    if M2 <: CompositeModel
        index_models!(index, right_run, M2)
    end
    # but we add right then left
    # need to do it like this to match `recursive_model_parse`
    !(M2 <: CompositeModel) && push!(index, right_run)
    !(M1 <: CompositeModel) && push!(index, left_run)
end

function index_models(model::Type{<:CompositeModel})
    index = Expr[]
    index_models!(index, :model, model)
    index
end

function index_models(::Type{<:AbstractSpectralModel})
    [:model]
end

function all_parameter_symbols(model::Type{<:AbstractSpectralModel})
    fieldnames(model)
end

function free_parameter_symbols(M::Type{<:AbstractSpectralModel})
    first(M.parameters[end].parameters)
end

function frozen_parameter_symbols(M::Type{<:AbstractSpectralModel})
    free_params = free_parameter_symbols(M)
    frozen = filter(i -> i ∉ free_params, all_parameter_symbols(M))
    tuple(frozen...)
end

function closure_parameter_symbols(::Type{<:AbstractSpectralModel})
    error("This specialisation should never need to be invoked.")
end

model_base_name(M::Type{<:AbstractSpectralModel}) = Base.typename(M).name

function _vector_to_named_tuple(params, names)
    statements = [
        :(params[$(i)]) for i in 1:length(names)
    ]
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

function all_parameters_to_named_tuple(params, model)
    names = all_parameter_symbols(model)
    _vector_to_named_tuple(params, names)
end

function all_parameters_to_named_tuple(model)
    names = all_parameter_symbols(model)
    statements = [
        :(getproperty(model, $(Meta.quot(s)))) for s in names
    ]
    :(NamedTuple{$(names)}(($(statements...),)))
end