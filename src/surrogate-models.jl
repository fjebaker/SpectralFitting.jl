
struct SurrogateSpectralModel{K,S,P,Z} <: AbstractSpectralModel
    surrogate::S
    params::P
    params_symbols::Z
    SurrogateSpectralModel(
        ::K,
        surrogate::S,
        params::P,
        params_symbols::Z,
    ) where {K,S,P,Z} = new{K,S,P,Z}(surrogate, params, params_symbols)
end

closurekind(::Type{<:SurrogateSpectralModel}) = WithClosures()
model_base_name(::Type{<:SurrogateSpectralModel{K}}) where {K} = :(SurrogateSpectralModel{$K})

# model generation
get_closure_param_fields(::Type{<:SurrogateSpectralModel}) = (:surrogate,)
get_param_types(::Type{<:SurrogateSpectralModel{K,S,P,Z}}) where {K,S,P,Z} = P.types
get_param_symbols(M::Type{<:SurrogateSpectralModel}) =
    [Symbol(:P, i) for i in eachindex(get_param_types(M))]

# runtime access
get_param_symbols(m::SurrogateSpectralModel) = m.params_symbols
function get_param(m::SurrogateSpectralModel, s::Symbol)
    i = findfirst(==(s), m.params_symbols)
    if !isnothing(i)
        m.params[i]
    else
        error("No such symbol: $s")
    end
end

modelkind(::Type{<:SurrogateSpectralModel{Additive}}) = Additive()
modelkind(::Type{<:SurrogateSpectralModel{Multiplicative}}) = Multiplicative()
modelkind(::Type{<:SurrogateSpectralModel{Convolutional}}) = Convolutional()

@fastmath function invoke!(
    flux,
    energy,
    ::Type{<:SurrogateSpectralModel{Multiplicative}},
    surrogate,
    params...,
)
    @inbounds for i in eachindex(flux)
        E = energy[i]
        v = (E, params...)
        flux[i] = surrogate(v)
    end
end

@fastmath function invoke!(
    flux,
    energy,
    ::Type{<:SurrogateSpectralModel{Additive}},
    surrogate,
    params...,
)
    integrate_over_flux!(flux, energy) do E
        surrogate((E, params...))
    end
end
