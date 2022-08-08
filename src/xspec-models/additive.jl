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

# broken
# symbol lookup error: lib/libXS.so: undefined symbol: _gfortran_string_len_trim
# @xspecmodel Additive :C_kerrdisk struct XS_KerrDisk{F1,F2,F3,F4,F5,F6,F7,F8,F9,F10}
#     "Normalisation."
#     K::F1 = FitParam(1.0)
#     "Rest frame line energy (keV)."
#     lineE::F2 = FrozenFitParam(6.4)
#     "Emissivity index for inner disk."
#     index1::F3 = FrozenFitParam(3.0)
#     "Emissivity index for outer disk."
#     index2::F4 = FrozenFitParam(3.0)
#     "Break radius seperating inner and outer disk (gᵣ)."
#     break_r::F5 = FrozenFitParam(6.0)
#     "Dimensionless black hole spin."
#     a::F6 = FitParam(0.998)
#     "Disk inclination angle to line of sight (degrees)."
#     incl::F7 = FrozenFitParam(30.0)
#     "Inner radius of the disk in units of rₘₛ."
#     inner_r::F8 = FrozenFitParam(1.0)
#     "Outer radius of the disk in units of rₘₛ."
#     outer_r::F9 = FrozenFitParam(400.0)
#     "Redshift."
#     z::F10 = FrozenFitParam(1.0)
# end

export XS_PowerLaw, XS_BlackBody, XS_BremsStrahlung
