@xspecmodel Convolutional :C_cflux struct XS_CalculateFlux{F1,F2,F3}
    "Minimum energy."
    E_min::F1 = FrozenFitParam(0.2)
    "Maximum energy."
    E_max::F2 = FrozenFitParam(2.0)
    "log (base 10) flux in erg / cm^2 / s"
    lg10Flux::F3 = FitParam(-10.0, lower_limit = -Inf, upper_limit = 0.0)
end

export XS_CalculateFlux
