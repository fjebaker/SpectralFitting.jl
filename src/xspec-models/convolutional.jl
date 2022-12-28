"""
    XS_CalculateFlux(E_min, E_max, lg10Flux)

$(FIELDS)
"""
@xspecmodel :C_cflux struct XS_CalculateFlux{T,F} <: AbstractSpectralModel{Convolutional}
    "Minimum energy."
    E_min::FitParam{T}
    "Maximum energy."
    E_max::FitParam{T}
    "log (base 10) flux in erg / cm^2 / s"
    log10Flux::FitParam{T}
    function XS_CalculateFlux(
        ; E_min = FitParam(0.2), E_max = FitParam(2.0), log10Flux = FitParam(-10.0, lower_limit = -Inf, upper_limit = 0.0)
    )
        new{parameter_type(E_min), FreeParameters{(:log10Flux,)}}(E_min, E_max, log10Flux)
    end
end

export XS_CalculateFlux
