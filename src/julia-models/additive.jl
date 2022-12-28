"""
    XS_PowerLaw(K, a)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, PowerLaw())
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
struct PowerLaw{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Photon index."
    a::FitParam{T}
    function PowerLaw(; K = FitParam(1.0), a = FitParam(2.0))
        new{parameter_type(K),FreeParameters{(:K, :a)}}(K, a)
    end
end
@inline function invoke!(flux, energy, ::Type{<:PowerLaw}, a)
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
invokemodel(energy, BlackBody())
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
struct BlackBody{T,F} <: AbstractSpectralModel{Additive}
    "Normalisation."
    K::FitParam{T}
    "Temperature (keV)."
    kT::FitParam{T}
    function BlackBody(; K = FitParam(1.0), kT = FitParam(3.0))
        new{parameter_type(K),FreeParameters{(:K, :kT)}}(K, kT)
    end
end
@inline function invoke!(flux, energy, ::Type{<:BlackBody}, kT)
    integration_kernel!(flux, energy) do E, δE
        8.0525 * E^2 * δE / (kT^4 * (exp(E / kT) - 1))
    end
end

export PowerLaw, BlackBody
