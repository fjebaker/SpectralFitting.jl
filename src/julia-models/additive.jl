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
struct PowerLaw{T,F} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Photon index."
    a::T
end
function PowerLaw(; K = FitParam(1.0), a = FitParam(2.0))
    PowerLaw{typeof(K),FreeParameters{(:K, :a)}}(K, a)
end
@inline @fastmath function invoke!(flux, energy, model::PowerLaw)
    let a = model.a
        α = 1 - a
        if (α ≈ 0.0)
            finite_diff_kernel!(flux, energy) do E
                log(E)
            end
        else
            α⁻¹ = inv(α)
            finite_diff_kernel!(flux, energy) do E
                α⁻¹ * E^α
            end
        end
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
struct BlackBody{T,F} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Temperature (keV)."
    kT::T
end
function BlackBody(; K = FitParam(1.0), kT = FitParam(3.0))
    BlackBody{typeof(K),FreeParameters{(:K, :kT)}}(K, kT)
end
@inline function invoke!(flux, energy, model::BlackBody)
    let kT = model.kT
        integration_kernel!(flux, energy) do E, δE
            8.0525 * E^2 * δE / (kT^4 * (exp(E / kT) - 1))
        end
    end
end

export PowerLaw, BlackBody
