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
                            AbstractSpectralModel{T,Multiplicative}
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::T
end
function XS_PhotoelectricAbsorption(; ηH = FitParam(1.0))
    XS_PhotoelectricAbsorption{typeof(ηH),FreeParameters{(:ηH,)}}(ηH)
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
@xspecmodel :C_wndabs struct XS_WarmAbsorption{T,F} <: AbstractSpectralModel{T,Multiplicative}
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::T
    "Window energy (keV)."
    Ew::T
end
function XS_WarmAbsorption(;
    ηH = FitParam(1.0),
    Ew = FitParam(1.0),
)
    XS_WarmAbsorption{typeof(ηH),FreeParameters{(:ηH, :Ew)}}(ηH, Ew)
end

export XS_PhotoelectricAbsorption, XS_WarmAbsorption
