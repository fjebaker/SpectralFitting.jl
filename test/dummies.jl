
struct DummyAdditive{T,F} <: AbstractSpectralModel{Additive}
    K::FitParam{T}
    a::FitParam{T}
    b::FitParam{T}
    function DummyAdditive(; K = FitParam(1.0), a = FitParam(1.0), b = FitParam(5.0))
        new{SpectralFitting.parameter_type(K),SpectralFitting.FreeParameters{(:K, :a)}}(
            K,
            a,
            b,
        )
    end
end
function SpectralFitting.invoke!(flux, energy, ::Type{<:DummyAdditive}, a, b)
    flux[:] .= a + b
end


struct DummyMultiplicative{T,F} <: AbstractSpectralModel{Multiplicative}
    a::FitParam{T}
    b::FitParam{T}
    function DummyMultiplicative(; a = FitParam(1.0), b = FitParam(5.0))
        new{SpectralFitting.parameter_type(a),SpectralFitting.FreeParameters{(:a,)}}(a, b)
    end
end
function SpectralFitting.invoke!(flux, energy, ::Type{<:DummyMultiplicative}, a, b)
    @. flux = flux * a + b
end
