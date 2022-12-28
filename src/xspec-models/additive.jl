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
@xspecmodel :C_powerlaw struct XS_PowerLaw{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Photon index."
    a::FitParam{T}
    function XS_PowerLaw(; K = FitParam(1.0), a = FitParam(1.0))
        new{SpectralFitting.parameter_type(K),SpectralFitting.FreeParameters{(:K, :a)}}(
            K,
            a,
        )
    end
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
@xspecmodel :C_bbody struct XS_BlackBody{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Temperature (keV)."
    T::FitParam{T}
    function XS_BlackBody(; K = FitParam(1.0), T = FitParam(3.0))
        new{parameter_type(K),FreeParameters{(:K, :T)}}(K, T)
    end
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
@xspecmodel :C_bremss struct XS_BremsStrahlung{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Plasma temperature (keV)."
    T::FitParam{T}
    function XS_BremsStrahlung(; K = FitParam(1.0), T = FitParam(7.0))
        new{parameter_type(K),FreeParameters{(:K, :T)}}(K, T)
    end
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
@xspecmodel :C_kerrdisk struct XS_KerrDisk{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Rest frame line energy (keV)."
    lineE::FitParam{T}
    "Emissivity index for inner disk."
    index1::FitParam{T}
    "Emissivity index for outer disk."
    index2::FitParam{T}
    "Break radius seperating inner and outer disk (gᵣ)."
    break_r::FitParam{T}
    "Dimensionless black hole spin."
    a::FitParam{T}
    "Disk inclination angle to line of sight (degrees)."
    incl::FitParam{T}
    "Inner radius of the disk in units of rₘₛ."
    inner_r::FitParam{T}
    "Outer radius of the disk in units of rₘₛ."
    outer_r::FitParam{T}
    "Redshift."
    z::FitParam{T}
    function XS_KerrDisk(;
        K = FitParam(1.0),
        lineE = FitParam(6.4),
        index1 = FitParam(3.0),
        index2 = FitParam(3.0),
        break_r = FitParam(6.0),
        a = FitParam(0.998, lower_limit = 0, upper_limit = 1.0),
        θ = FitParam(30.0),
        inner_r = FitParam(1.0),
        outer_r = FitParam(400.0),
        z = FitParam(0.0),
    )
        new{parameter_type(K),FreeParameters{(:K, :a)}}(
            K,
            lineE,
            index1,
            index2,
            break_r,
            a,
            θ,
            inner_r,
            outer_r,
            z,
        )
    end
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
@xspecmodel :C_kyrline struct XS_KyrLine{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Dimensionless black hole spin."
    a::FitParam{T}
    "Observer inclination (0 is on pole, degrees)."
    θ_obs::FitParam{T}
    "Inner radius of the disk in units of GM/c²"
    inner_r::FitParam{T}
    "0: integrate from rᵢₙ. 1: integrate from rₘₛ."
    ms_flag::FitParam{T}
    "Outer radius of the disk in units of GM/c²"
    outer_r::FitParam{T}
    "Rest frame line energy (keV)."
    lineE::FitParam{T}
    α::FitParam{T}
    β::FitParam{T}
    "Break radius seperating inner and outer disk (GM/c²)."
    break_r::FitParam{T}
    "Overall Doppler shift."
    z::FitParam{T}
    "0: isotropic emission, 1: Laor's limb darkening, 2: Haard's limb brightening."
    limb::FitParam{T}
    function XS_KerrDisk(;
        K = FitParam(1.0),
        a = FitParam(0.998, lower_limit = 0, upper_limit = 1.0),
        θ = FitParam(30.0),
        inner_r = FitParam(1.0),
        outer_r = FitParam(400.0),
        lineE = FitParam(6.4),
        α = FitParam(3.0),
        β = FitParam(3.0),
        break_r = FitParam(6.0),
        z = FitParam(0.0),
        limb = FitParam(1.0),
    )
        new{parameter_type(K),FreeParameters{(:K, :a, :θ)}}(
            K,
            a,
            θ,
            inner_r,
            outer_r,
            lineE,
            α,
            β,
            break_r,
            z,
            limb,
        )
    end
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
@xspecmodel :C_laor struct XS_Laor{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Rest frame line energy (keV)."
    lineE::FitParam{T}
    "Power law dependence of emissitivy. Scales R⁻ᵅ."
    a::FitParam{T}
    "Inner radius of the accretion disk (GM/c)."
    inner_r::FitParam{T}
    "Outer radius of the accretion disk (GM/c)."
    outer_r::FitParam{T}
    "Disk inclination angle to line of sight (degrees, 0 is pole on)."
    θ::FitParam{T}
    function XS_Laor(;
        K = FitParam(1.0),
        lineE = FitParam(6.4),
        a = FitParam(3.0),
        inner_r = FitParam(1.235),
        outer_r = FitParam(400.0),
        θ = FitParam(30.0, upper_limit = 180),
    )
        new{parameter_type(K),FreeParameters{(:K, :lineE)}}(
            K,
            lineE,
            a,
            inner_r,
            outer_r,
            θ,
        )
    end
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
@xspecmodel :C_diskline struct XS_DiskLine{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Rest frame line energy (keV)."
    lineE::FitParam{T}
    "Power law dependence of emissitivy. If < 10, scales Rᵅ."
    β::FitParam{T}
    "Inner radius of the accretion disk (GM/c)."
    inner_r::FitParam{T}
    "Outer radius of the accretion disk (GM/c)."
    outer_r::FitParam{T}
    "Disk inclination angle to line of sight (degrees, 0 is pole on)."
    θ::FitParam{T}
    function XS_DiskLine(;
        K = FitParam(1.0),
        lineE = FitParam(6.7),
        β = FitParam(-2.0),
        inner_r = FitParam(10.0),
        outer_r = FitParam(1000.0),
        θ = FitParam(30.0, upper_limit = 180),
    )
        new{parameter_type(K),FreeParameters{(:K, :lineE)}}(
            K,
            lineE,
            β,
            inner_r,
            outer_r,
            θ,
        )
    end
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
@xspecmodel :C_gaussian struct XS_Gaussian{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation"
    K::FitParam{T}
    "Line wavelength in Angstrom."
    E::FitParam{T}
    "Line width in Angstrom."
    σ::FitParam{T}
    function XS_Gaussian(; K = FitParam(1.0), E = FitParam(6.4), σ = FitParam(1.0))
        new{parameter_type(K),FreeParameters{(:K, :σ)}}(K, E, σ)
    end
end

export XS_PowerLaw,
    XS_BlackBody,
    XS_BremsStrahlung,
    XS_Laor,
    XS_DiskLine,
    XS_KerrDisk,
    XS_KyrLine,
    XS_Gaussian
