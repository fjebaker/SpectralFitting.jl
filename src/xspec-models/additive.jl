@xspecmodel Additive :C_powerlaw struct XS_PowerLaw{F}
    "Normalisation."
    K::F = FitParameter(1.0)
    "Photon index."
    a::F = FitParameter(0.5)
end


export XS_PowerLaw
