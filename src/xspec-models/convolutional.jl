"""
    XS_CalculateFlux(E_min, E_max, lg10Flux)

$(FIELDS)
"""
@xspecmodel :C_cflux struct XS_CalculateFlux{T} <: AbstractSpectralModel{T,Convolutional}
    "Minimum energy."
    E_min::T
    "Maximum energy."
    E_max::T
    "log (base 10) flux in erg / cm^2 / s"
    log10Flux::T
end
function XS_CalculateFlux(;
    E_min = FitParam(0.2, frozen = true),
    E_max = FitParam(2.0, frozen = true),
    log10Flux = FitParam(-10.0, lower_limit = -Inf, upper_limit = 0.0),
)
    XS_CalculateFlux(E_min, E_max, log10Flux)
end

export XS_CalculateFlux
