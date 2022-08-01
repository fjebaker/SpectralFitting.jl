# @xspecmodel :C_thcomp struct XS_ThermalCompton{Convolutional}
#     "Normalisation."
#     K::Float64 = 1.0
#     "Asymptotic power-law photon index."
#     Γ::Float64 = 1.0
# end

@xspecmodel :C_cflux struct XS_ConvFlux{Convolutional}
    "Minimum energy."
    E_min::Float64 = 0.2
    "Maximum energy."
    E_max::Float64 = 2.0
    "log (base 10) flux in erg / cm^2 / s"
    lg10Flux::Float64 = -10.0
end

export XS_ConvFlux
