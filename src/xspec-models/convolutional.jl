@xspecmodel Convolutional :C_cflux struct XS_CalculateFlux{F}
    "Minimum energy."
    E_min::F = FitParameter(0.2, frozen = true)
    "Maximum energy."
    E_max::F = FitParameter(2.0, frozen = true)
    "log (base 10) flux in erg / cm^2 / s"
    lg10Flux::F = FitParameter(-10.0, lower_bound = -Inf, upper_bound = 0.0)
end

export XS_CalculateFlux
