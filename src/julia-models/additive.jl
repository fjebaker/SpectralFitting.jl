@with_kw struct PowerLaw{F} <: AbstractSpectralModel
    K::F = FitParam(1.0)
    a::F = FitParam(0.5)
end
modelkind(::Type{<:PowerLaw}) = Additive()
@fastmath function invoke!(flux, energy, ::Type{<:PowerLaw}, a)
    α = 1 - a
    α⁻¹ = inv(α)
    integrate_over_flux!(flux, energy) do E
        α⁻¹ * E^α
    end
end

@with_kw struct BlackBody{F} <: AbstractSpectralModel
    K::F = FitParam(1.0)
    kT::F = FitParam(2.0)
end
modelkind(::Type{<:BlackBody}) = Additive()
@fastmath function invoke!(flux, energy, ::Type{<:BlackBody}, kT)
    integrate_over_flux!(flux, energy) do E
        8.0525 * E^2 / (kT^4 * (exp(E / kT) - 1))
    end
end

export PowerLaw
