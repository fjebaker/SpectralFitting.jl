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
@xspecmodel Multiplicative :C_phabs struct XS_PhotoelectricAbsorption{F}
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::F = FitParam(1.0)
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
@xspecmodel Multiplicative :C_wndabs struct XS_WarmAbsorption{F1,F2}
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::F1 = FitParam(1.0)
    "Window energy (keV)."
    Ew::F2 = FitParam(1.0)
end

export XS_PhotoelectricAbsorption, XS_WarmAbsorption
