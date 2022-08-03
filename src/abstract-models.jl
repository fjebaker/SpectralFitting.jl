export AbstractSpectralModel, modelkind, invokemodel!, value, upperbound, lowerbound, as_distribution

# traits
abstract type AbstractSpectralModelKind end
struct Multiplicative <: AbstractSpectralModelKind end
struct Additive <: AbstractSpectralModelKind end
struct Convolutional <: AbstractSpectralModelKind end

# fitting parameters
abstract type AbstractFitParameter end
value(f::AbstractFitParameter) = f.val
value(x::Number) = x
setvalue!(f::AbstractFitParameter, val) = f.val = val
upperbound(f::AbstractFitParameter) = f.upper_bound
lowerbound(f::AbstractFitParameter) = f.lower_bound
setlowerbound!(f::AbstractFitParameter, lb) = f.lower_bound = lb
setupperbound!(f::AbstractFitParameter, ub) = f.upper_bound = ub
isfrozen(f::AbstractFitParameter) = f.frozen
freeze!(f::AbstractFitParameter) = f.frozen = true

unfreeze!(f::AbstractFitParameter) = f.frozen = false
as_distribution(f::AbstractFitParameter) = error("Not implemented for $(typeof(f)).")

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
        flux ./= value(m.K)
    else
        invoke!(flux, energy, m)
    end
end
