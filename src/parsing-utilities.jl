unpack_model(::Type{<:SpectralFitting.CompositeSpectralModel{M1,M2,O}}) where {M1,M2,O} =
    (M1, M2, O)
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
