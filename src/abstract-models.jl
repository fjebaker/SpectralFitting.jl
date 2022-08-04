export AbstractSpectralModel,
    modelkind,
    invokemodel!,
    value,
    upperbound,
    lowerbound,
    as_distribution,
    AbstractSpectralModelKind,
    Multiplicative,
    Additive,
    Convolutional

# traits
abstract type AbstractSpectralModelKind end
struct Multiplicative <: AbstractSpectralModelKind end
struct Additive <: AbstractSpectralModelKind end
struct Convolutional <: AbstractSpectralModelKind end

abstract type AbstractModelImplementation end
struct XSPECImplementation <: AbstractModelImplementation end
struct JuliaImplementation <: AbstractModelImplementation end

abstract type AbstractModelconstructionkind end
struct NonTrivialConstruction <: AbstractModelconstructionkind end
struct TrivialConstruction <: AbstractModelconstructionkind end

# models
abstract type AbstractSpectralModel end
invoke!(flux, energy, m::AbstractSpectralModel) = error("Not defined for $(typeof(m)).")
modelkind(m::Type{<:AbstractSpectralModel}) = error("Not defined for $(typeof(m)).")
implementation(::Type{<:AbstractSpectralModel}) = JuliaImplementation()
constructionkind(::Type{<:AbstractSpectralModel}) = TrivialConstruction()

function modelinfo(m::M) where {M<:AbstractSpectralModel}
    params = join([get_value(p) for p in get_all_parameters(m)], ", ")
    "$(model_type_name(M))[$(params)]"
end

numbertype(::AbstractSpectralModel) = Sys.WORD_SIZE == 64 ? Float64 : Float32
# needed:
get_parameter_symbols(::M) where {M<:AbstractSpectralModel} = fieldnames(M)
get_parameter(m::AbstractSpectralModel, s::Symbol) = getproperty(m, s)
# optional:
get_all_parameters(m::M) where {M<:AbstractSpectralModel} =
    (get_parameter(m, p) for p in get_parameter_symbols(m))
# todo: make this a proper iterator? also better name
parameter_symbol_pairs(m::M) where {M<:AbstractSpectralModel} =
    (p => get_parameter(m, p) for p in get_parameter_symbols(m))
model_type_name(::Type{M}) where {M<:AbstractSpectralModel} = Base.typename(M).name

invokemodel!(flux, energy, m::M) where {M<:AbstractSpectralModel} =
    invokemodel!(flux, energy, m, modelkind(M))
function invokemodel!(flux, energy, model, ::K) where {K}
    invoke!(flux, energy, model)
    flux
end
@fastmath function invokemodel!(flux, energy, model, ::Additive)
    invoke!(flux, energy, model)
    flux ./= get_value(model.K)
end

function Base.show(io::IO, ::MIME"text/plain", m::M) where {M<:AbstractSpectralModel}
    params = [String(s) => p for (s, p) in parameter_symbol_pairs(m)]
    print(io, "$(model_type_name(M))\n")

    pad = maximum(i -> length(first(i)), params) + 1

    for (s, val) in params
        print(io, "   $(rpad(s, pad)) => ")
        println(io, val)
    end
end
