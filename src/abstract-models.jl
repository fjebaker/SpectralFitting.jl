export AbstractSpectralModel, modelkind, invokemodel!

# traits
abstract type AbstractSpectralModelKind end
struct Multiplicative <: AbstractSpectralModelKind end
struct Additive <: AbstractSpectralModelKind end
struct Convolutional <: AbstractSpectralModelKind end

# fitting paramters
abstract type AbstractFitParameter end
value(f::AbstractFitParameter) = f.val
value(x::Number) = x
upperbound(f::AbstractFitParameter) = f.ub
lowerbound(f::AbstractFitParameter) = f.lb
isfrozen(f::AbstractFitParameter) = f.frozen
freeze!(f::AbstractFitParameter) = f.frozen = true
unfreeze!(f::AbstractFitParameter) = f.frozen = false

# models
abstract type AbstractSpectralModel end
modelkind(m::Type{<:AbstractSpectralModel}) = error("Not defined for $(typeof(m)).")
function modelinfo(m::M) where {M<:AbstractSpectralModel}
    params = join([value(getproperty(m, p)) for p in fieldnames(M)], ", ")
    "$(Base.typename(M).name)[$(params)]"
end
invoke!(flux, energy, m::AbstractSpectralModel) = error("Not defined for $(typeof(m)).")
# bootstrap this so we can do things like additive normalisation
# this could also be generated
function invokemodel!(flux, energy, m::M) where {M<:AbstractSpectralModel}
    if modelkind(M) == Additive
        invoke!(flux, energy, m)
        flux .*= value(m.K)
    else
        invoke!(flux, energy, m)
    end
end
