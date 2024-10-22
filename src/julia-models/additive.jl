"""
    PowerLaw

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
struct PowerLaw{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Photon index."
    a::T
end
function PowerLaw(; K = FitParam(1.0), a = FitParam(2.0))
    PowerLaw(K, a)
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
    BlackBody

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
struct BlackBody{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Temperature (keV)."
    kT::T
end
function BlackBody(; K = FitParam(1.0), kT = FitParam(3.0))
    BlackBody{typeof(K)}(K, kT)
end
@inline function invoke!(flux, energy, model::BlackBody)
    let kT = model.kT
        integration_kernel!(flux, energy) do E, δE
            8.0525 * E^2 * δE / (kT^4 * (exp(E / kT) - 1))
        end
    end
end

"""
    GaussianLine

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, GuassianLine())
```

```
                       GaussianLine                
        ┌────────────────────────────────────────┐ 
   0.09 │                                        │ 
        │            ..                          │ 
        │           .':                          │ 
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
        │        :       '.                      │ 
      0 │.......:         :......................│ 
        └────────────────────────────────────────┘ 
         0                                     20  
                          E (keV)                  
```
"""
struct GaussianLine{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Mean"
    μ::T
    "Standard deviation"
    σ::T
end
function GaussianLine(; K = FitParam(1.0), μ = FitParam(6.4), σ = FitParam(1.0))
    GaussianLine{typeof(K)}(K, μ, σ)
end
@inline function invoke!(flux, energy, model::GaussianLine)
    let μ = model.μ, σ = model.σ
        integration_kernel!(flux, energy) do Elow, δE
            E = Elow + δE / 2
            δE * inv(σ * √(2π)) * exp(-1 * (E - μ)^2 / (2 * σ^2))
        end
    end
end

"""
    DeltaLine

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, DeltaLine())
```

```
                        DeltaLine                 
       ┌────────────────────────────────────────┐ 
   0.4 │         .                              │ 
       │         :                              │ 
       │         :                              │ 
       │         :                              │ 
       │         :                              │ 
       │         :                              │ 
       │         ::                             │ 
       │         ::                             │ 
       │         ::                             │ 
       │         ::                             │ 
       │         ::                             │ 
       │         ::                             │ 
       │         ::                             │ 
       │         ::                             │ 
     0 │.........::.............................│ 
       └────────────────────────────────────────┘ 
        0                                     20  
                         E (keV)                 
```

!!! note
    The `DeltaLine` model is not a true delta function, as this would be
    extremely difficult to define in a numerical model that needs to be able to
    propagate gradients. Instead, it is a very narrow [`GaussianLine`](@ref)
    model.
"""
struct DeltaLine{W<:Number,T} <: AbstractSpectralModel{T,Additive}
    _width::W
    "Normalisation."
    K::T
    "Energy at which the delta function spikes."
    E::T
end

function DeltaLine(; K = FitParam(1.0), E = FitParam(5.0), width = 1e-2)
    DeltaLine{typeof(width),typeof(K)}(width, K, E)
end

Reflection.get_closure_symbols(::Type{<:DeltaLine}) = (:_width,)

Reflection.get_parameter_symbols(model::Type{<:DeltaLine}) = fieldnames(model)[2:end]

@inline function invoke!(flux, energy, model::DeltaLine{T}) where {T}
    # we can't actually have a diract delta because that would ruin
    # the ability to run dual numbers through the system. What we can do instead
    # is have a miniscule gaussian
    invoke!(flux, energy, GaussianLine(promote(model.K, model.E, model._width)...))
end

export PowerLaw, BlackBody, GaussianLine, DeltaLine
