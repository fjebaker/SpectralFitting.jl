
# standard julia models for testing

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

# table models for testing

struct DummyAdditiveTableModel{T,D,F} <: AbstractTableModel{D,Additive}
    table::D
    K::FitParam{T}
    a::FitParam{T}
    b::FitParam{T}
    function DummyAdditiveTableModel(
        K::FitParam{T},
        a::FitParam{T},
        b::FitParam{T},
    ) where {T}
        # table is just an interpolation anywhere so lambda for tests
        table = (x) -> x^2
        new{T,typeof(table),SpectralFitting.FreeParameters{(:K, :a)}}(table, K, a, b)
    end
end
function SpectralFitting.invoke!(
    flux,
    energy,
    ::Type{<:DummyAdditiveTableModel},
    a,
    b,
    table,
)
    flux[:] .= table(a) + b
end
function DummyAdditiveTableModel(; K = FitParam(1.0), a = FitParam(1.0), b = FitParam(2.0))
    DummyAdditiveTableModel(K, a, b)
end

struct DummyMultiplicativeTableModel{T,D,F} <: AbstractTableModel{D,Multiplicative}
    table::D
    a::FitParam{T}
    b::FitParam{T}
    function DummyMultiplicativeTableModel(a::FitParam{T}, b::FitParam{T}) where {T}
        # table is just an interpolation anywhere so lambda for tests
        table = (x, k) -> k * x
        new{T,typeof(table),SpectralFitting.FreeParameters{(:a,)}}(table, a, b)
    end
end
function SpectralFitting.invoke!(
    flux,
    energy,
    ::Type{<:DummyMultiplicativeTableModel},
    a,
    b,
    table,
)
    @. flux = table(flux, a) + b
end
function DummyMultiplicativeTableModel(; a = FitParam(1.0), b = FitParam(2.0))
    DummyMultiplicativeTableModel(a, b)
end
