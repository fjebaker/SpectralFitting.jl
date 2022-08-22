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
