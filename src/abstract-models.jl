export AbstractSpectralModel, modelkind, invokemodel!, value, upperbound, lowerbound, as_distribution

# traits
abstract type AbstractSpectralModelKind end
struct Multiplicative <: AbstractSpectralModelKind end
struct Additive <: AbstractSpectralModelKind end
struct Convolutional <: AbstractSpectralModelKind end

abstract type AbstractModelImplementation end
struct XSPECImplementation <: AbstractModelImplementation end
struct JuliaImplementation <: AbstractModelImplementation end

# models
abstract type AbstractSpectralModel end
invoke!(flux, energy, m::AbstractSpectralModel) = error("Not defined for $(typeof(m)).")
modelkind(m::Type{<:AbstractSpectralModel}) = error("Not defined for $(typeof(m)).")
implementation(m::Type{<:AbstractSpectralModel}) = JuliaImplementation()

function modelinfo(m::M) where {M<:AbstractSpectralModel}
    params = join([get_value(getproperty(m, p)) for p in fieldnames(M)], ", ")
    "$(Base.typename(M).name)[$(params)]"
end

numbertype(m::AbstractSpectralModel) = Sys.WORD_SIZE == 64 ? Float64 : Float32

invokemodel!(flux, energy, m::M) where {M<:AbstractSpectralModel} = invokemodel!(flux, energy, m, modelkind(M))
function invokemodel!(flux, energy, model, ::K) where {K}
    invoke!(flux, energy, model)
    flux
end
@fastmath function invokemodel!(flux, energy, model, ::Additive)
    invoke!(flux, energy, model)
    flux ./= get_value(model.K)
end

function Base.show(io::IO, ::MIME"text/plain", m::M) where {M<:AbstractSpectralModel}
    params = [String(p) => getproperty(m, p) for p in fieldnames(M)]
    print(io, "$(Base.typename(M).name)\n")

    pad = maximum(i -> length(first(i)), params) + 1

    for (s, val) in params
        print(io, "   $(rpad(s, pad)) => ")
        println(io, val)
    end
end
