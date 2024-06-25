"""
    abstract type AbstractModelWrapper{M,T,K} <: AbstractSpectralModel{T,K} end

First field of the struct must be `model`.
"""
abstract type AbstractModelWrapper{M<:AbstractSpectralModel,T,K} <:
              AbstractSpectralModel{T,K} end

normalisation(model::AbstractModelWrapper{M,T,Additive}) where {M,T} =
    normalisation(model.model)

function Reflection.get_closure_symbols(
    M::Type{<:AbstractModelWrapper{Model}},
) where {Model}
    # we ignore the `model` field, since that will be given by the constructor
    (fieldnames(M)[2:end]..., Reflection.get_closure_symbols(Model)...)
end

Reflection.get_parameter_symbols(::Type{<:AbstractModelWrapper{M}}) where {M} =
    Reflection.get_parameter_symbols(M)

function Reflection.closure_parameter_lenses(
    M::Type{<:AbstractModelWrapper},
    info::Reflection.ModelInfo,
)
    num_closures = fieldcount(M) - 1 # ignore the `model` field

    my_closures = map(info.closure_symbols[1:num_closures]) do s
        :(getfield($(info.lens), $(Meta.quot(s))))
    end
    model_closures = map(info.closure_symbols[num_closures+1:end]) do s
        :(getfield($(info.lens).model, $(Meta.quot(s))))
    end
    vcat(my_closures, model_closures)
end

function Reflection.parameter_lenses(
    ::Type{<:AbstractModelWrapper},
    info::Reflection.ModelInfo,
)
    map(info.symbols) do s
        :(getfield($(info.lens).model, $(Meta.quot(s))))
    end
end

function Reflection.make_constructor(
    M::Type{<:AbstractModelWrapper{Model}},
    closures::Vector,
    params::Vector,
    T::Type,
) where {Model}
    num_closures = fieldcount(M) - 1 # ignore the `model` field
    my_closures = closures[1:num_closures]
    model_constructor =
        Reflection.make_constructor(Model, closures[num_closures+1:end], params, T)
    :($(Base.typename(M).name)($(model_constructor), $(my_closures...)))
end
