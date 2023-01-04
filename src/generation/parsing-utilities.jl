struct ModelInfo
    symbols::Vector{Symbol}
    free::Vector{Symbol}
    frozen::Vector{Symbol}
    generated_symbols::Vector{Symbol}
    lens::Union{Symbol,Expr}
    type::Type
end

function getinfo(model::Type{<:AbstractSpectralModel}; lens::Union{Symbol,Expr}=:(model))
    symbs = [all_parameter_symbols(model)...]
    free = [free_parameter_symbols(model)...]
    frozen = [frozen_parameter_symbols(model)...]
    ModelInfo(symbs, free, frozen, [Base.gensym(s) for s in symbs], lens, model)
end

function getinfo(model::Type{<:CompositeModel})
    infos = ModelInfo[]
    _addinfo!(infos, model, :(model))
    infos
end

function _addinfo!(infos, model::Type{<:AbstractSpectralModel}, lens::Union{Symbol,Expr})
    push!(infos, getinfo(model; lens = lens))
end
function _addinfo!(infos, model::Type{<:CompositeModel}, lens::Union{Symbol,Expr})
    right_run = :(getproperty($lens, :right))
    _addinfo!(infos, model.parameters[2], right_run)
    left_run = :(getproperty($lens, :left))
    _addinfo!(infos, model.parameters[1], left_run)
end

unpack_model(::Type{<:CompositeModel{M1,M2,O}}) where {M1,M2,O} = (M1, M2, O)
unpack_model(m::CompositeModel{M1,M2,O}) where {M1,M2,O} =
    m.left, m.right, operation_symbol(O)

@inline function _recursive_model_parse(callback, model)
    left_model, right_model, operator = unpack_model(model)
    left = recursive_model_parse(callback, left_model)
    right = recursive_model_parse(callback, right_model)
    callback((left, right, operator))
end

recursive_model_parse(_, model::Type{<:AbstractSpectralModel}) = model
function recursive_model_parse(callback, model::Type{<:CompositeModel})
    _recursive_model_parse(callback, model)
    # Nothing
end

recursive_model_parse(_, model::AbstractSpectralModel) = model
function recursive_model_parse(callback, model::CompositeModel)
    _recursive_model_parse(callback, model)
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

function all_parameters_to_named_tuple(params, model)
    names = all_parameter_symbols(model)
    _vector_to_named_tuple(params, names)
end

function all_parameters_to_named_tuple(model)
    names = all_parameter_symbols(model)
    statements = [:(getproperty(model, $(Meta.quot(s)))) for s in names]
    :(NamedTuple{$(names)}(($(statements...),)))
end
