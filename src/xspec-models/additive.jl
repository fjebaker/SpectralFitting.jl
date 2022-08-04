@xspecmodel Additive :C_powerlaw struct XS_PowerLaw{F}
    "Normalisation."
    K::F = FitParam(1.0)
    "Photon index."
    a::F = FitParam(0.5)
end

@xspecmodel Additive :C_bbody struct XS_BlackBody{F}
    "Normalisation."
    K::F = FitParam(1.0)
    "Temperature (keV)"
    T::F = FitParam(3.0)
end

@xspecmodel Additive :C_bremss struct XS_BremsStrahlung{F}
    "Normalisation."
    K::F = FitParam(1.0)
    "Plasma temperature (keV)"
    T::F = FitParam(7.0)
end

export XS_PowerLaw, XS_BlackBody, XS_BremsStrahlung
