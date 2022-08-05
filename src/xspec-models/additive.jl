@xspecmodel Additive :C_powerlaw struct XS_PowerLaw{F1,F2}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Photon index."
    a::F2 = FitParam(0.5)
end

@xspecmodel Additive :C_bbody struct XS_BlackBody{F1,F2}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Temperature (keV)"
    T::F2 = FitParam(3.0)
end

@xspecmodel Additive :C_bremss struct XS_BremsStrahlung{F1,F2}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Plasma temperature (keV)"
    T::F2 = FitParam(7.0)
end

export XS_PowerLaw, XS_BlackBody, XS_BremsStrahlung
