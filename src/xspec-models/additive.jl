@xspecmodel :C_compbb struct XS_ComptonizationBlackBody{Additive}
    "Normalisation."
    K::Float64 = 1.0
    "Blackbody temperature (keV)."
    kT::Float64 = 1.0
    "Electron temperature of hot plasma (keV)."
    kTe::Float64 = 50.0
    "Optical depth of the plasma."
    τ::Float64 = 0.1
end

@xspecmodel :C_bremss struct XS_ThermalBrems{Additive}
    "Normalisation."
    K::Float64 = 1.0
    "Plasma temperature (keV)."
    pT::Float64 = 0.5
end

@xspecmodel :C_powerlaw struct XS_PowerLaw{Additive}
    "Normalisation."
    K::Float64 = 1.0
    "Photon index."
    a::Float64 = 0.5
end

@xspecmodel :C_zpowerlw struct XS_RedshiftPowerLaw{Additive}
    "Normalisation."
    K::Float64 = 1.0
    "Photon index."
    a::Float64 = 0.5
    "Redshift."
    z::Float64 = 1.0
end

@xspecmodel :C_srcut struct XS_SynchrotronCutoff{Additive}
    "Normalisation."
    K::Float64 = 1.0
    "Radio spectral index."
    a::Float64 = 0.5
    "Break Hz: approx frequency at which flux has dropped by factor of 10 from powerlaw."
    cutoff::Float64 = 2.42e17
end

export XS_ComptonizationBlackBody,
    XS_ThermalBrems, XS_PowerLaw, XS_RedshiftPowerLaw, XS_SynchrotronCutoff
