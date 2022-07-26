abstract type AbstractModelKind end

struct Multiplicative <: AbstractModelKind end
struct Additive <: AbstractModelKind end
struct Convolutional <: AbstractModelKind end

abstract type AbstractSpectralModel{K<:AbstractModelKind} end

struct CompositeSpectralModel{K1,K2,M1,M2} <: AbstractSpectralModel{K2}
    m1::M1
    m2::M2
    function CompositeSpectralModel(m1::AbstractSpectralModel, m2::AbstractSpectralModel)
        new{modelkind(m1), modelkind(m2), typeof(m1), typeof(m2)}(m1, m2)
    end
end

function prettystring(::T) where {T<:AbstractSpectralModel{K}} where {K}
    Base.typename(T).name
end

# function prettystring(::T) where {T<:AbstractSpectralModel{Additive}}
#     Base.typename(T).name
# end

function prettystring(m::CompositeSpectralModel{Multiplicative,Additive})
    n1 = prettystring(m.m1)
    n2 = prettystring(m.m2)
    "[$n1 * $n2]"
end

function prettystring(m::CompositeSpectralModel{Multiplicative,Multiplicative})
    n1 = prettystring(m.m1)
    n2 = prettystring(m.m2)
    "$n1 * $n2"
end

function prettystring(m::CompositeSpectralModel{Additive,Additive})
    n1 = prettystring(m.m1)
    n2 = prettystring(m.m2)
    "$n1 + $n2"
end

function Base.show(io::IO, m::CompositeSpectralModel{K1,K2}) where {K1,K2}
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == CompositeSpectralModel
        print(io, "CompositeModel($(modelkind(m)))")
    else
        model_string = prettystring(m)
        print(io, "CompositeSpectralModel ($K2):\n  $model_string")
    end
end

# interface
function invokemodel(input, m::AbstractSpectralModel; kwargs...)
    output = zeros(eltype(input), size(input))
    err = zeros(eltype(input), size(input))
    invokemodel!(output, err, input, m; kwargs...)
    output
end

function invokemodel!(output, input, m::AbstractSpectralModel; kwargs...)
    err = zeros(eltype(input), size(input))
    invokemodel!(output, err, input, m; kwargs...)
    output
end

@inline invokemodel!(output, err, input, m::AbstractSpectralModel; kwargs...) =
    error("Not defined for model $m.")

# apply model
(m::AbstractSpectralModel)(energy::AbstractArray; kwargs...) =
    invokemodel(energy, m; kwargs...)
(m::AbstractSpectralModel{Additive})(energy::AbstractArray; kwargs...) =
    m.K .* invokemodel(energy, m; kwargs...)
(m::CompositeSpectralModel{Multiplicative,K})(energy::AbstractArray; kwargs...) where {K} =
    m.m1(energy; kwargs...) .* m.m2(energy; kwargs...)
(m::CompositeSpectralModel{Additive,Additive})(energy::AbstractArray; kwargs...) =
    m.m1(energy; kwargs...) .+ m.m2(energy; kwargs...)


# don't think this interface should be supported, in favour of using the above
# to trial out if a model looks like it should, and then favour a "build model" command,
# which creates optimised code
# # apply model with alternative parameters
# (m::AbstractSpectralModel)(energy::AbstractArray, params; kwargs...) =
#     invokemodel(typeof(m)(params...), energy; kwargs...)
# function (m::AbstractSpectralModel{Additive})(energy::AbstractArray, params; kwargs...)
#     m2 = typeof(m)(params...)
#     m2.K .* invokemodel(m2, energy; kwargs...)
# end

# get parameters of the model
modelparams(m::AbstractSpectralModel) = [getproperty(m, fn) for fn in fieldnames(typeof(m))]

# get normalisation of Additive models
knorm(::AbstractSpectralModel) = error("Only defined for Additive models.")
knorm(m::AbstractSpectralModel{Additive}) = m.K

modelkind(cm::CompositeSpectralModel{K1, K2}) where {K1, K2} = modelkind(cm.m2)
modelkind(::AbstractSpectralModel{K}) where {K} = K

# model algebra
Base.:*(::AbstractSpectralModel, ::AbstractSpectralModel) = error("Left model must be Multiplicative")
Base.:*(m1::AbstractSpectralModel{Multiplicative}, m2::AbstractSpectralModel) = CompositeSpectralModel(m1, m2)

Base.:+(::AbstractSpectralModel, ::AbstractSpectralModel) = error("Both models must be Additive.")
Base.:+(m1::AbstractSpectralModel{Additive}, m2::AbstractSpectralModel{Additive}) = CompositeSpectralModel(m1, m2)

# Base.:∘()

export AbstractModelKind,
    AbstractSpectralModel,
    Multiplicative,
    Additive,
    Convolutional,
    params,
    modelkind
