export AbstractSpectralModel,
    AbstractSpectralModelKind,
    Multiplicative,
    Additive,
    Convolutional,
    modelkind,
    AbstractSpectralModelImplementation,
    XSPECImplementation,
    JuliaImplementation,
    implementation,
    AbstractSpectralModelClosureType,
    WithClosures,
    WithoutClosures,
    closurekind,
    has_closure_params,
    #get_param_symbols,
    #get_param_types,
    #model_base_name,
    #get_closure_param_fields,
    get_param_symbols,
    get_param,
    get_param_count,
    get_all_params,
    get_all_params_by_value,
    get_param_symbol_pairs,
    invokemodel!,
    flux_count


# models
abstract type AbstractSpectralModel end
numbertype(::AbstractSpectralModel) = Sys.WORD_SIZE == 64 ? Float64 : Float32

# traits
abstract type AbstractSpectralModelKind end
struct Multiplicative <: AbstractSpectralModelKind end
struct Additive <: AbstractSpectralModelKind end
struct Convolutional <: AbstractSpectralModelKind end

modelkind(m::Type{<:AbstractSpectralModel}) = error("Not defined for $(typeof(m)).")

abstract type AbstractSpectralModelImplementation end
struct XSPECImplementation <: AbstractSpectralModelImplementation end
struct JuliaImplementation <: AbstractSpectralModelImplementation end

implementation(::Type{<:AbstractSpectralModel}) = JuliaImplementation()

abstract type AbstractSpectralModelClosureType end
struct WithClosures <: AbstractSpectralModelClosureType end
struct WithoutClosures <: AbstractSpectralModelClosureType end

closurekind(::Type{<:AbstractSpectralModel}) = WithoutClosures()

has_closure_params(::WithClosures) = true
has_closure_params(::WithoutClosures) = false
has_closure_params(M::Type{<:AbstractSpectralModel}) = has_closure_params(closurekind(M))
has_closure_params(::M) where {M<:AbstractSpectralModel} = has_closure_params(M)

# implementation interface
# never to be called directly
# favour `invokemodel!` instead
invoke!(flux, energy, ::Type{M}, params...) where {M<:AbstractSpectralModel} =
    error("Not defined for $(M).")

# optional interface
get_param_symbols(M::Type{<:AbstractSpectralModel}) = fieldnames(M)
get_param_types(M::Type{<:AbstractSpectralModel}) = M.types

model_base_name(M::Type{<:AbstractSpectralModel}) = Base.typename(M).name

# only needed for WithClosures()
get_closure_param_fields(::Type{<:AbstractSpectralModel}) = ()

# minimal param accessors
get_param_types(::M) where {M<:AbstractSpectralModel} = get_param_types(M)
get_param_symbols(::M) where {M<:AbstractSpectralModel} = get_param_symbols(M)
get_param(m::AbstractSpectralModel, s::Symbol) = getproperty(m, s)

# emergent param accessors
get_param_count(M::Type{<:AbstractSpectralModel}) = length(get_param_types(M))

get_all_params(m::AbstractSpectralModel) where {M} =
    (get_param(m, p) for p in get_param_symbols(m))
get_all_params_by_value(m::AbstractSpectralModel) =
    (get_value(i) for i in get_all_params(m))
# todo: make this a proper iterator? also better name
get_param_symbol_pairs(m::M) where {M<:AbstractSpectralModel} =
    (p => get_param(m, p) for p in get_param_symbols(m))

# invokation wrappers
invokemodel!(f, e, m::M) where {M<:AbstractSpectralModel} =
    invokemodel!(f, e, M, get_all_params_by_value(m)...)
invokemodel!(f, e, ::Type{M}, p...) where {M<:AbstractSpectralModel} =
    invokemodel!(f, e, modelkind(M), M, p...)
@fastmath function invokemodel!(
    flux,
    energy,
    ::Additive,
    M::Type{<:AbstractSpectralModel},
    K,
    p...,
)
    invoke!(flux, energy, M, p...)
    flux ./= K
end
@fastmath function invokemodel!(
    flux,
    energy,
    ::AbstractSpectralModelKind,
    M::Type{<:AbstractSpectralModel},
    p...,
)
    invoke!(flux, energy, M, p...)
    flux
end

# bindings to generated functions

flux_count(model::AbstractSpectralModel) = generated_maximum_flux_count(model)

# printing

function modelinfo(m::M) where {M<:AbstractSpectralModel}
    params = join([get_value(p) for p in get_all_params(m)], ", ")
    "$(model_base_name(M))[$(params)]"
end

function Base.show(io::IO, ::MIME"text/plain", m::M) where {M<:AbstractSpectralModel}
    params = [String(s) => p for (s, p) in get_param_symbol_pairs(m)]
    print(io, "$(model_base_name(M))\n")

    pad = maximum(i -> length(first(i)), params) + 1

    for (s, val) in params
        print(io, "   $(rpad(s, pad)) => ")
        println(io, val)
    end
end
