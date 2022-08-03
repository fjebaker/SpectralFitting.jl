export AbstractSpectralModel, modelkind, invokemodel!, value, upperbound, lowerbound, as_distribution

# traits
abstract type AbstractSpectralModelKind end
struct Multiplicative <: AbstractSpectralModelKind end
struct Additive <: AbstractSpectralModelKind end
struct Convolutional <: AbstractSpectralModelKind end

# models
abstract type AbstractSpectralModel end
modelkind(m::Type{<:AbstractSpectralModel}) = error("Not defined for $(typeof(m)).")
function modelinfo(m::M) where {M<:AbstractSpectralModel}
    params = join([get_value(getproperty(m, p)) for p in fieldnames(M)], ", ")
    "$(Base.typename(M).name)[$(params)]"
end
numbertype(m::AbstractSpectralModel) = Sys.WORD_SIZE == 64 ? Float64 : Float32
invoke!(flux, energy, m::AbstractSpectralModel) = error("Not defined for $(typeof(m)).")
# bootstrap this so we can do things like additive normalisation
# this could also be generated
function invokemodel!(flux, energy, m::M) where {M<:AbstractSpectralModel}
    if modelkind(M) == Additive
        invoke!(flux, energy, m)
        flux ./= get_value(m.K)
    else
        invoke!(flux, energy, m)
    end
end
