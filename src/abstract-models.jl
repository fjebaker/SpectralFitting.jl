abstract type AbstractModelKind end

struct Multiplicative <: AbstractModelKind end
struct Additive <: AbstractModelKind end
struct Convolutional <: AbstractModelKind end
struct Composite{K1,K2} <: AbstractModelKind
    mk1::K1
    mk2::K2
end

abstract type AbstractSpectralModel{K<:AbstractModelKind} end

# interface
invokemodel(m::AbstractSpectralModel, ::AbstractArray; kwargs...) =
    error("Not defined for model $m.")

# apply model
(m::AbstractSpectralModel)(energy::AbstractArray; kwargs...) =
    invokemodel(m, energy; kwargs...)
(m::AbstractSpectralModel{Additive})(energy::AbstractArray; kwargs...) =
    m.K .* invokemodel(m, energy; kwargs...)

# apply model with alternative parameters
(m::AbstractSpectralModel)(energy::AbstractArray, params; kwargs...) =
    invokemodel(typeof(m)(params), energy; kwargs...)
function (m::AbstractSpectralModel{Additive})(energy::AbstractArray, params; kwargs...)
    m2 = typeof(m)(params)
    m2.K .* invokemodel(m2, energy; kwargs...)
end

# get parameters of the model
modelparams(m::AbstractSpectralModel) = [getproperty(m, fn) for fn in fieldnames(typeof(m))]

# get normalisation of Additive models
knorm(::AbstractSpectralModel) = error("Only defined for Additive models.")
knorm(m::AbstractSpectralModel{Additive}) = m.K

# get right operation kind
# M(f)
rmodelkind(K::AbstractModelKind) = K
rmodelkind(::Composite{K1,K2}) where {K1,K2} = K2
rmodelkind(::AbstractSpectralModel{K}) where {K} = rmodelkind(K)

# get left operation kind
# f(M)
lmodelkind(K::AbstractModelKind) = K
lmodelkind(::Composite{K1,K2}) where {K1,K2} = K1
lmodelkind(::AbstractSpectralModel{K}) where {K} = rmodelkind(K)

modelkind(::AbstractSpectralModel{K}) where {K} = K


export AbstractModelKind,
    AbstractSpectralModel,
    Multiplicative,
    Additive,
    Convolutional,
    params,
    rmodelkind,
    lmodelkind
