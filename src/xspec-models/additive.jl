"""
    XS_PowerLaw(K, a)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_PowerLaw())
```

```
                      XS_PowerLaw
       ┌────────────────────────────────────────┐
   0.5 │                                        │
       │:                                       │
       │:                                       │
       │:                                       │
       │:                                       │
       │:                                       │
       │:                                       │
       │ :                                      │
       │ :                                      │
       │  :                                     │
       │   :.                                   │
       │    ':..                                │
       │        ''':......                      │
       │                  ''''''''''''''........│
     0 │                                        │
       └────────────────────────────────────────┘
        0                                     20
                         E (keV)
```
"""
@xspecmodel Additive :C_powerlaw struct XS_PowerLaw{F1,F2}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Photon index."
    a::F2 = FitParam(0.5)
end

"""
    XS_BlackBody(K, T)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_BlackBody())
```

```
                      XS_BlackBody
       ┌────────────────────────────────────────┐
   0.2 │                                        │
       │                                        │
       │                                        │
       │                                        │
       │                                        │
       │                                        │
       │      .:''':..                          │
       │    .:       ''.                        │
       │   .'           ':.                     │
       │   :              ''..                  │
       │  :                  ':.                │
       │ :                     '':.             │
       │.:                         ''..         │
       │:                              '':....  │
     0 │'                                    '''│
       └────────────────────────────────────────┘
        0                                     20
                         E (keV)
```
"""
@xspecmodel Additive :C_bbody struct XS_BlackBody{F1,F2}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Temperature (keV)."
    T::F2 = FitParam(3.0)
end

"""
    XS_BremsStrahlung(K, T)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_BremsStrahlung())
```

```
                  XS_BremsStrahlung
     ┌────────────────────────────────────────┐
   2 │                                        │
     │.                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │'.                                      │
     │ :                                      │
   0 │  ':....................................│
     └────────────────────────────────────────┘
      0                                     20
                       E (keV)
```
"""
@xspecmodel Additive :C_bremss struct XS_BremsStrahlung{F1,F2}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Plasma temperature (keV)."
    T::F2 = FitParam(7.0)
end


"""
    XS_KerrDisk(K, lineE, index1, index2, break_r, a, incl, inner_r, outer_r)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_KerrDisk())
```

```
                        XS_KerrDisk
        ┌────────────────────────────────────────┐
   0.05 │                                        │
        │                                        │
        │                                .       │
        │                               .:       │
        │                              :::       │
        │                            .:' '.      │
        │                           .:    :      │
        │                         ..'     :      │
        │                         :'      :      │
        │                       .'        :      │
        │                     .:'         :      │
        │                  .:''           :      │
        │               .::'              :      │
        │            ..:'                 :      │
      0 │.........:'''                    :......│
        └────────────────────────────────────────┘
         0                                      8
                          E (keV)
```
"""
@xspecmodel Additive :C_kerrdisk struct XS_KerrDisk{F1,F2,F3,F4,F5,F6,F7,F8,F9,F10}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Rest frame line energy (keV)."
    lineE::F2 = FrozenFitParam(6.4)
    "Emissivity index for inner disk."
    index1::F3 = FrozenFitParam(3.0)
    "Emissivity index for outer disk."
    index2::F4 = FrozenFitParam(3.0)
    "Break radius seperating inner and outer disk (gᵣ)."
    break_r::F5 = FrozenFitParam(6.0)
    "Dimensionless black hole spin."
    a::F6 = FitParam(0.998, upper_limit = 1.0)
    "Disk inclination angle to line of sight (degrees)."
    incl::F7 = FrozenFitParam(30.0)
    "Inner radius of the disk in units of rₘₛ."
    inner_r::F8 = FrozenFitParam(1.0)
    "Outer radius of the disk in units of rₘₛ."
    outer_r::F9 = FrozenFitParam(400.0)
    "Redshift."
    z::F10 = FrozenFitParam(0.0)
end
register_model_data(XS_KerrDisk, "kerrtable.fits")


"""
    XS_KyrLine(K, a, θ_obs, inner_r, ms_flag, outer_r, lineE, α, β, break_r, z, limb)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_KyrLine())
```

```
                        XS_KyrLine
        ┌────────────────────────────────────────┐
   0.05 │                                        │
        │                                        │
        │                                :       │
        │                                :.      │
        │                              :.':      │
        │                             :'  :      │
        │                            :    :      │
        │                          .'     :      │
        │                        .:'      :      │
        │                       .'        :      │
        │                     .:          :      │
        │                   .'            :      │
        │                .:'              :      │
        │            ..:'                 :      │
      0 │.........:'''                    :......│
        └────────────────────────────────────────┘
         0                                      8
                          E (keV)
```
"""
@xspecmodel Additive :C_kyrline struct XS_KyrLine{F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Dimensionless black hole spin."
    a::F2 = FitParam(0.998, upper_limit = 1.0)
    "Observer inclination (0 is on pole, degrees)."
    θ_obs::F3 = FitParam(30.0)
    "Inner radius of the disk in units of GM/c²"
    inner_r::F4 = FrozenFitParam(1.0)
    "0: integrate from rᵢₙ. 1: integrate from rₘₛ."
    ms_flag::F5 = FrozenFitParam(1)
    "Outer radius of the disk in units of GM/c²"
    outer_r::F6 = FrozenFitParam(400.0)
    "Rest frame line energy (keV)."
    lineE::F7 = FrozenFitParam(6.4)
    α::F8 = FrozenFitParam(3.0)
    β::F9 = FrozenFitParam(3.0)
    "Break radius seperating inner and outer disk (GM/c²)."
    break_r::F10 = FrozenFitParam(400.0)
    "Overall Doppler shift."
    z::F11 = FrozenFitParam(0.0)
    "0: isotropic emission, 1: Laor's limb darkening, 2: Haard's limb brightening."
    limb::F12 = FrozenFitParam(1)
end
register_model_data(XS_KyrLine, "KBHline01.fits")


"""
    XS_Laor(K, lineE, a, inner_r, outer_r, incl)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 10.0, 100))
invokemodel(energy, XS_Laor())
```

```
                          XS_Laor
        ┌────────────────────────────────────────┐
   0.06 │                                        │
        │                                        │
        │                         ::             │
        │                         ::             │
        │                        : :             │
        │                       :  :             │
        │                      :   :             │
        │                     :'   :             │
        │                   .:     :             │
        │                  :'      :             │
        │                .'        :             │
        │              .:'         :             │
        │            ..'           :             │
        │          .:'              :            │
      0 │.......:''                 :............│
        └────────────────────────────────────────┘
         0                                     10
                          E (keV)
```
"""
@xspecmodel Additive :C_laor struct XS_Laor{F1,F2,F3,F4,F5,F6}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Rest frame line energy (keV)."
    lineE::F2 = FitParam(6.4)
    "Power law dependence of emissitivy. Scales R⁻ᵅ."
    a::F3 = FrozenFitParam(3.0)
    "Inner radius of the accretion disk (GM/c)."
    inner_r::F4 = FrozenFitParam(1.235)
    "Outer radius of the accretion disk (GM/c)."
    outer_r::F5 = FrozenFitParam(400.0)
    "Disk inclination angle to line of sight (degrees, 0 is pole on)."
    incl::F6 = FitParam(30.0, upper_limit = 180)
end
register_model_data(XS_Laor, "ari.mod")

"""
    XS_DiskLine(K, lineE, β, inner_r, outer_r, incl)

$(FIELDS)

# Example

```julia
energy = collect(range(4.0, 8.0, 100))
invokemodel(energy, XS_DiskLine())
```

```
                        XS_DiskLine
        ┌────────────────────────────────────────┐
   0.09 │                                        │
        │                           .            │
        │                           :            │
        │                           ::           │
        │                         . ::           │
        │                         : ::           │
        │                         :'':           │
        │                        .'  :           │
        │                        :    :          │
        │                        :    :          │
        │                       .'    :          │
        │                       :     :          │
        │                     .:      '.         │
        │                   .:'        :         │
      0 │...............:'''           :.........│
        └────────────────────────────────────────┘
         4                                      8
                          E (keV)
```
"""
@xspecmodel Additive :C_diskline struct XS_DiskLine{F1,F2,F3,F4,F5,F6}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Rest frame line energy (keV)."
    lineE::F2 = FitParam(6.7)
    "Power law dependence of emissitivy. If < 10, scales Rᵅ."
    β::F3 = FrozenFitParam(-2.0)
    "Inner radius of the accretion disk (GM/c)."
    inner_r::F4 = FrozenFitParam(10.0)
    "Outer radius of the accretion disk (GM/c)."
    outer_r::F5 = FrozenFitParam(1000.0)
    "Disk inclination angle to line of sight (degrees, 0 is pole on)."
    incl::F6 = FitParam(30.0, upper_limit = 180)
end

"""
    XS_Gaussian(K, E, σ)

$(FIELDS)

# Example

```julia
energy = collect(range(4.0, 8.0, 100))
invokemodel(energy, XS_Gaussian())
```

```
                        XS_Gaussian                
        ┌────────────────────────────────────────┐ 
   0.09 │                                        │ 
        │            .                           │ 
        │           : :                          │ 
        │           : :                          │ 
        │           : '.                         │ 
        │          .'  :                         │ 
        │          :   :                         │ 
        │          :   :                         │ 
        │          :   '.                        │ 
        │         :     :                        │ 
        │         :     :                        │ 
        │         :     :                        │ 
        │        .'      :                       │ 
        │        :       :                       │ 
      0 │.......:         :......................│ 
        └────────────────────────────────────────┘ 
         0                                     20  
                          E (keV)                  
```
"""
@xspecmodel Additive :C_gaussian struct XS_Gaussian{F1,F2,F3}
    "Normalisation"
    K::F1 = FitParam(1.0)
    "Line wavelength in Angstrom."
    E::F2 = FrozenFitParam(6.4)
    "Line width in Angstrom."
    σ::F3 = FitParam(1.0)
end

export XS_PowerLaw,
    XS_BlackBody,
    XS_BremsStrahlung,
    XS_Laor,
    XS_DiskLine,
    XS_KerrDisk,
    XS_KyrLine,
    XS_Gaussian
