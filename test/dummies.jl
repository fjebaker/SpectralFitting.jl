
# standard julia models for testing
struct DummyAdditive{T,F} <: AbstractSpectralModel{T,Additive}
    K::T
    a::T
    b::T
end
function DummyAdditive(; K = FitParam(1.0), a = FitParam(1.0), b = FitParam(5.0))
    DummyAdditive{typeof(K),SpectralFitting.FreeParameters{(:K, :a)}}(K, a, b)
end
function SpectralFitting.invoke!(flux, energy, model::DummyAdditive)
    let a = model.a, b = model.b
        flux[:] .= a + b
    end
end

struct DummyMultiplicative{T,F} <: AbstractSpectralModel{T,Multiplicative}
    a::T
    b::T
end
function DummyMultiplicative(; a = FitParam(1.0), b = FitParam(5.0))
    DummyMultiplicative{typeof(a),SpectralFitting.FreeParameters{(:a,)}}(a, b)
end
function SpectralFitting.invoke!(flux, energy, model::DummyMultiplicative)
    let a = model.a, b = model.b
        @. flux = flux * a + b
    end
end

# table models for testing

struct DummyAdditiveTableModel{D,T,F} <: AbstractTableModel{T,Additive}
    table::D
    K::T
    a::T
    b::T
end
function DummyAdditiveTableModel(
    K::T,
    a::T,
    b::T
) where {T}
    # table is just an interpolation anywhere so lambda for tests
    table = (x) -> x^2
    DummyAdditiveTableModel{typeof(table),T,SpectralFitting.FreeParameters{(:K, :a)}}(table, K, a, b)
end
function SpectralFitting.invoke!(
    flux,
    energy,
    model::DummyAdditiveTableModel
)
    let a = model.a, b = model.b, table = model.table
        flux[:] .= table(a) + b
    end
end
function DummyAdditiveTableModel(; K = FitParam(1.0), a = FitParam(1.0), b = FitParam(2.0))
    DummyAdditiveTableModel(K, a, b)
end

struct DummyMultiplicativeTableModel{D,T,F} <: AbstractTableModel{T,Multiplicative}
    table::D
    a::T
    b::T
end
function DummyMultiplicativeTableModel(a::T, b::T) where {T}
    # table is just an interpolation anywhere so lambda for tests
    table = (x, k) -> k * x
    DummyMultiplicativeTableModel{typeof(table),T,SpectralFitting.FreeParameters{(:a,)}}(table, a, b)
end
function SpectralFitting.invoke!(
    flux,
    energy,
    model::DummyMultiplicativeTableModel
)
    let a = model.a, b = model.b, table = model.table
        @. flux = table(flux, a) + b
    end
end
function DummyMultiplicativeTableModel(; a = FitParam(1.0), b = FitParam(2.0))
    DummyMultiplicativeTableModel(a, b)
end
