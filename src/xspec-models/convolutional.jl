@xspecmodel Convolutional :C_cflux struct XS_CalculateFlux{F}
    "Minimum energy."
    E_min::F = FitParam(0.2, frozen = true, error = 0.0)
    "Maximum energy."
    E_max::F = FitParam(2.0, frozen = true, error = 0.0)
    "log (base 10) flux in erg / cm^2 / s"
    lg10Flux::F = FitParam(-10.0, lower_limit = -Inf, upper_limit = 0.0)
end

export XS_CalculateFlux
