"""
    XS_PowerLaw(K, a)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokeflux(energy, PowerLaw())
```

```
                        PowerLaw
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
@with_kw struct PowerLaw{F1,F2} <: AbstractSpectralModel
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Photon index."
    a::F2 = FitParam(0.5)
end
modelkind(::Type{<:PowerLaw}) = Additive()
@fastmath function invoke!(flux, energy, ::Type{<:PowerLaw}, a)
    α = 1 - a
    α⁻¹ = inv(α)
    finite_diff_kernel!(flux, energy) do E
        α⁻¹ * E^α
    end
end

"""
    XS_BlackBody(K, T)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokeflux(energy, BlackBody())
```

```
                        BlackBody
       ┌────────────────────────────────────────┐
   0.2 │                                        │
       │                                        │
       │                                        │
       │                                        │
       │                                        │
       │                                        │
       │      .:''':..                          │
       │     :'      '':.                       │
       │   .'           ':.                     │
       │  .:               '..                  │
       │  :                  ':.                │
       │ .'                     ':..            │
       │ :                         ''...        │
       │:                              '''....  │
     0 │:                                    '''│
       └────────────────────────────────────────┘
        0                                     20
                         E (keV)
```
"""
@with_kw struct BlackBody{F1,F2} <: AbstractSpectralModel
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Temperature (keV)."
    kT::F2 = FitParam(3.0)
end
modelkind(::Type{<:BlackBody}) = Additive()
@fastmath function invoke!(flux, energy, ::Type{<:BlackBody}, kT)
    integration_kernel!(flux, energy) do E, δE
        8.0525 * E^2 * δE / (kT^4 * (exp(E / kT) - 1))
    end
end

export PowerLaw, BlackBody
