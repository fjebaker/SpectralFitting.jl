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
@xspecmodel :C_powerlaw struct XS_PowerLaw{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Photon index."
    a::T
end
function XS_PowerLaw(; K = FitParam(1.0), a = FitParam(1.0))
    XS_PowerLaw(K, a)
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
@xspecmodel :C_bbody struct XS_BlackBody{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Temperature (keV)."
    T::T
end
function XS_BlackBody(; K = FitParam(1.0), T = FitParam(3.0))
    XS_BlackBody(K, T)
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
@xspecmodel :C_bremss struct XS_BremsStrahlung{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Plasma temperature (keV)."
    T::T
end
function XS_BremsStrahlung(; K = FitParam(1.0), T = FitParam(7.0))
    XS_BremsStrahlung(K, T)
end


"""
    XS_KerrDisk(K, lineE, index1, index2, break_r, a, θ, inner_r, outer_r)

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
@xspecmodel :C_kerrdisk struct XS_KerrDisk{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Rest frame line energy (keV)."
    lineE::T
    "Emissivity index for inner disk."
    index1::T
    "Emissivity index for outer disk."
    index2::T
    "Break radius seperating inner and outer disk (gᵣ)."
    break_r::T
    "Dimensionless black hole spin."
    a::T
    "Disk inclination angle to line of sight (degrees)."
    θ::T
    "Inner radius of the disk in units of rₘₛ."
    inner_r::T
    "Outer radius of the disk in units of rₘₛ."
    outer_r::T
    "Redshift."
    z::T
end
function XS_KerrDisk(;
    K = FitParam(1.0),
    lineE = FitParam(6.4, frozen = true),
    index1 = FitParam(3.0, frozen = true),
    index2 = FitParam(3.0, frozen = true),
    break_r = FitParam(6.0, frozen = true),
    a = FitParam(0.998, lower_limit = 0, upper_limit = 0.998),
    θ = FitParam(30.0, lower_limit = 0, upper_limit = 90.0),
    inner_r = FitParam(1.0, frozen = true),
    outer_r = FitParam(400.0, frozen = true),
    z = FitParam(0.0, frozen = true),
)
    XS_KerrDisk(K, lineE, index1, index2, break_r, a, θ, inner_r, outer_r, z)
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
@xspecmodel :C_kyrline struct XS_KyrLine{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Dimensionless black hole spin."
    a::T
    "Observer inclination (0 is on pole, degrees)."
    θ::T
    "Inner radius of the disk in units of GM/c²"
    inner_r::T
    "0: integrate from rᵢₙ. 1: integrate from rₘₛ."
    ms_flag::T
    "Outer radius of the disk in units of GM/c²"
    outer_r::T
    "Rest frame line energy (keV)."
    lineE::T
    α::T
    β::T
    "Break radius seperating inner and outer disk (GM/c²)."
    break_r::T
    "Overall Doppler shift."
    z::T
    "0: isotropic emission, 1: Laor's limb darkening, 2: Haard's limb brightening."
    limb::T
end
function XS_KyrLine(;
    K = FitParam(1.0),
    a = FitParam(0.998, lower_limit = 0, upper_limit = 1.0),
    θ = FitParam(30.0),
    inner_r = FitParam(1.0, frozen = true),
    ms_flag = FitParam(0.0, frozen = true),
    outer_r = FitParam(400.0, frozen = true),
    lineE = FitParam(6.4, frozen = true),
    α = FitParam(3.0, frozen = true),
    β = FitParam(3.0, frozen = true),
    break_r = FitParam(6.0, frozen = true),
    z = FitParam(0.0, frozen = true),
    limb = FitParam(1.0, frozen = true),
)
    XS_KyrLine(K, a, θ, inner_r, ms_flag, outer_r, lineE, α, β, break_r, z, limb)
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
@xspecmodel :C_laor struct XS_Laor{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Rest frame line energy (keV)."
    lineE::T
    "Power law dependence of emissitivy. Scales R⁻ᵅ."
    a::T
    "Inner radius of the accretion disk (GM/c)."
    inner_r::T
    "Outer radius of the accretion disk (GM/c)."
    outer_r::T
    "Disk inclination angle to line of sight (degrees, 0 is pole on)."
    θ::T
end
function XS_Laor(;
    K = FitParam(1.0),
    lineE = FitParam(6.4),
    a = FitParam(3.0, frozen = true),
    inner_r = FitParam(1.235, frozen = true),
    outer_r = FitParam(400.0, frozen = true),
    θ = FitParam(30.0, upper_limit = 180, frozen = true),
)
    XS_Laor(K, lineE, a, inner_r, outer_r, θ)
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
@xspecmodel :C_diskline struct XS_DiskLine{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Rest frame line energy (keV)."
    lineE::T
    "Power law dependence of emissitivy. If < 10, scales Rᵅ."
    β::T
    "Inner radius of the accretion disk (GM/c)."
    inner_r::T
    "Outer radius of the accretion disk (GM/c)."
    outer_r::T
    "Disk inclination angle to line of sight (degrees, 0 is pole on)."
    θ::T
end
function XS_DiskLine(;
    K = FitParam(1.0),
    lineE = FitParam(6.7),
    β = FitParam(-2.0, frozen = true),
    inner_r = FitParam(10.0, frozen = true),
    outer_r = FitParam(1000.0, frozen = true),
    θ = FitParam(30.0, upper_limit = 180, frozen = true),
)
    XS_DiskLine(K, lineE, β, inner_r, outer_r, θ)
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
@xspecmodel :C_gaussian struct XS_Gaussian{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
    "Line wavelength in Angstrom."
    E::T
    "Line width in Angstrom."
    σ::T
end
function XS_Gaussian(;
    K = FitParam(1.0),
    E = FitParam(6.4, frozen = true),
    σ = FitParam(1.0),
)
    XS_Gaussian(K, E, σ)
end

export XS_PowerLaw,
    XS_BlackBody,
    XS_BremsStrahlung,
    XS_Laor,
    XS_DiskLine,
    XS_KerrDisk,
    XS_KyrLine,
    XS_Gaussian
