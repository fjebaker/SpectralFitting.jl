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

invoke!(flux, energy, ::Type{M}, params...) where {M<:AbstractSpectralModel} =
    error("Not defined for $(M).")
# invoke!(flux, energy, m::AbstractSpectralModel) = error("Not defined for $(typeof(m)).")
modelkind(m::Type{<:AbstractSpectralModel}) = error("Not defined for $(typeof(m)).")
implementation(::Type{<:AbstractSpectralModel}) = JuliaImplementation()
constructionkind(::Type{<:AbstractSpectralModel}) = TrivialConstruction()

is_trivially_constructed(::NonTrivialConstruction) = false
is_trivially_constructed(::TrivialConstruction) = true
is_trivially_constructed(::M) where {M<:AbstractSpectralModel} =
    is_trivially_constructed(constructionkind(M))

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

get_all_parameters_by_value(m::AbstractSpectralModel) =
    (get_value(i) for i in get_all_parameters(m))
# todo: make this a proper iterator? also better name
parameter_symbol_pairs(m::M) where {M<:AbstractSpectralModel} =
    (p => get_parameter(m, p) for p in get_parameter_symbols(m))
model_type_name(::Type{M}) where {M<:AbstractSpectralModel} = Base.typename(M).name

invokemodel!(f, e, m::M) where {M<:AbstractSpectralModel} =
    invokemodel!(f, e, M, get_all_parameters_by_value(m)...)
invokemodel!(f, e, ::Type{M}, p...) where {M<:AbstractSpectralModel} =
    invokemodel!(f, e, modelkind(M), M, p...)
@fastmath function invokemodel!(
    flux,
    energy,
    ::Additive,
    ::Type{M},
    K,
    p...,
) where {M<:AbstractSpectralModel}
    invoke!(flux, energy, M, p...)
    flux ./= K
end
@fastmath function invokemodel!(
    flux,
    energy,
    ::AbstractSpectralModelKind,
    ::Type{M},
    p...,
) where {M<:AbstractSpectralModel}
    invoke!(flux, energy, M, p...)
    flux
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
