unpack_model(::Type{<:CompositeSpectralModel{M1,M2,O}}) where {M1,M2,O} = (M1, M2, O)
unpack_model(m::CompositeSpectralModel{M1,M2,O}) where {M1,M2,O} =
    m.left, m.right, operation_symbol(O)

@inline function __recursive_model_parse(callback, model)
    left_model, right_model, operator = unpack_model(model)
    left = recursive_model_parse(callback, left_model)
    right = recursive_model_parse(callback, right_model)
    callback((left, right, operator))
end

recursive_model_parse(_, model::Type{<:AbstractSpectralModel}) = model
function recursive_model_parse(callback, model::Type{<:CompositeSpectralModel})
    __recursive_model_parse(callback, model)
    Nothing
end

recursive_model_parse(_, model::AbstractSpectralModel) = model
function recursive_model_parse(callback, model::CompositeSpectralModel)
    __recursive_model_parse(callback, model)
    nothing
end

function make_unique_readable_symbol(p, symbol_bank; delim = '_')
    i = 1
    symb = Symbol(p, delim, i)
    while symb in symbol_bank
        i += 1
        symb = Symbol(p, delim, i)
    end
    symb
end

function index_models!(
    index,
    running,
    ::Type{<:CompositeSpectralModel{M1,M2}},
) where {M1,M2}
    left_run = :(getproperty($running, :left))
    right_run = :(getproperty($running, :right))
    if M2 <: CompositeSpectralModel
        index_models!(index, right_run, M2)
    else
        push!(index, right_run)
    end
    if M1 <: CompositeSpectralModel
        index_models!(index, left_run, M1)
    else
        push!(index, left_run)
    end
end

function index_models(model::Type{<:CompositeSpectralModel})
    index = Expr[]
    index_models!(index, :model, model)
    index
end


function index_models(::Type{<:AbstractSpectralModel})
    [:model]
end
