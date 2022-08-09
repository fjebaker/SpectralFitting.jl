
@xspecmodel Multiplicative :C_phabs struct XS_PhotoelectricAbsorption{F}
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::F = FitParam(1.0)
end

@xspecmodel Multiplicative :C_wndabs struct XS_WarmAbsorption{F1,F2}
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::F1 = FitParam(1.0)
    "Window energy (keV)."
    Ew::F2 = FitParam(1.0)
end

export XS_PhotoelectricAbsorption, XS_WarmAbsorption
