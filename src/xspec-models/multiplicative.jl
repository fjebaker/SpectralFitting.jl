"""
    XS_PhotoelectricAbsorption(ηH)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_PhotoelectricAbsorption())
```

```
             XS_PhotoelectricAbsorption
     ┌────────────────────────────────────────┐
   1 │       ...''''''''''''''''''''''''''''''│
     │      .'                                │
     │     :                                  │
     │    :'                                  │
     │    :                                   │
     │   :                                    │
     │   :                                    │
     │   :                                    │
     │  :                                     │
     │  :                                     │
     │  :                                     │
     │  :                                     │
     │  :                                     │
     │ :                                      │
   0 │.:                                      │
     └────────────────────────────────────────┘
      0                                     20
                       E (keV)
```
"""
@xspecmodel :C_phabs struct XS_PhotoelectricAbsorption{T,F} <:
                            AbstractSpectralModel{Multiplicative}
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::FitParam{T}
    function XS_PhotoelectricAbsorption(; ηH = FitParam(1.0))
        new{parameter_type(ηH),FreeParameters{(:ηH,)}}(ηH)
    end
end

"""
    XS_WarmAbsorption(ηH, Ew)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_WarmAbsorption())
```

```
                    XS_WarmAbsorption
       ┌────────────────────────────────────────┐
     1 │':      ...''':'''''''''''''''''''''''''│
       │ :    .:'                               │
       │ :   .'                                 │
       │ :  .:                                  │
       │ :  :                                   │
       │ :  :                                   │
       │ : :                                    │
       │ : :                                    │
       │ : :                                    │
       │ : :                                    │
       │ ::                                     │
       │ ::                                     │
       │  :                                     │
       │  :                                     │
   0.2 │  :                                     │
       └────────────────────────────────────────┘
        0                                     20
                         E (keV)
```
"""
@xspecmodel :C_wndabs struct XS_WarmAbsorption{T,F} <: AbstractSpectralModel{Multiplicative}
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::FitParam{T}
    "Window energy (keV)."
    Ew::FitParam{T}
    function XS_WarmAbsorption(;
        ηH = FitParam(1.0),
        Ew = FitParam(1.0),
    )
        new{parameter_type(ηH),FreeParameters{(:ηH, :Ew)}}(ηH, Ew)
    end
end

export XS_PhotoelectricAbsorption, XS_WarmAbsorption
