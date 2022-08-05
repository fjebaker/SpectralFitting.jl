
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

additional_invoke_parameters(m::SurrogateSpectralModel) = (m.surrogate,)
make_additional_invoke_parameters_symbols(::Type{<:SurrogateSpectralModel}, symb) =
    (Symbol(symb, '_', :surrogate),)

make_additional_invoke_parameters_symbols(T::Type{<:AbstractSpectralModel}, symb) =
    error("Not implemented for $T.")
make_additional_invoke_parameters_symbols(m::M, symb) where {M<:AbstractSpectralModel} =
    make_additional_invoke_parameters_symbols(M, symb)

function assemble_invoke(flux, symb, M::Type{<:SurrogateSpectralModel}, params)
    kind = typeof(modelkind(M))
    surrogate_param_symbol = first(make_additional_invoke_parameters_symbols(M, symb))
    :(invokemodel!(
        $flux,
        energy,
        SurrogateSpectralModel{$kind},
        $surrogate_param_symbol,
        $(params...),
    ))
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
