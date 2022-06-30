
@with_kw struct PowerLaw <: AbstractSpectralModel{Additive}
    K::Float64 = 1.0
    a::Float64 = 0.5
end

function invokemodel(m::PowerLaw, energy::AbstractArray)
    energy .^ (-m.a)
end

export PowerLaw
