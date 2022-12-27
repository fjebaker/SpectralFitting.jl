using Test
using SpectralFitting

# test specific setup
struct DummyAdditive{T,F} <: AbstractSpectralModel
    K::FitParam{T}
    a::FitParam{T}
    b::FitParam{T}
    function DummyAdditive(;K = FitParam(1.0), a=FitParam(1.0), b=FitParam(5.0))
        new{SpectralFitting.parameter_type(K),SpectralFitting.FreeParameters{(:K,:a)}}(K, a, b)
    end
end
SpectralFitting.modelkind(::Type{<:DummyAdditive{T,F}}) where {T,F} = Additive()

add_model = DummyAdditive()

SpectralFitting.modelkind(add_model)
SpectralFitting.modelkind(typeof(add_model))

new_model = add_model + add_model
# SpectralFitting.generated_model_types(new_model)
# SpectralFitting.generated_model_types(PowerLaw() + PowerLaw())
# info = SpectralFitting.assemble_aggregate_info(typeof(new_model))
