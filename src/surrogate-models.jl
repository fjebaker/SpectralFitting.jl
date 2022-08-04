
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

constructionkind(::Type{<:SurrogateSpectralModel}) = NonTrivialConstruction()
get_parameter_symbols(m::SurrogateSpectralModel) = m.params_symbols
function get_parameter(m::SurrogateSpectralModel, s::Symbol)
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

@fastmath function invoke!(flux, energy, m::SurrogateSpectralModel)
    @inbounds for i in eachindex(flux)
        E = energy[i]
        v = (E, m.params...)
        flux[i] = m.surrogate(v)
    end
end
