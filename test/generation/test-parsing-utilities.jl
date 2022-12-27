using Test
using SpectralFitting

# test specific setup
struct DummyModel_TPU{T,F} <: AbstractSpectralModel
    K::FitParam{T}
    a::FitParam{T}
    b::FitParam{T}
    function DummyModel_TPU(;K = FitParam(1.0), a=FitParam(1.0), b=FitParam(5.0))
        new{SpectralFitting.parameter_type(K),SpectralFitting.FreeParameters{(:K,:a)}}(K, a, b)
    end
end
SpectralFitting.modelkind(::Type{<:DummyModel_TPU}) = Additive()

model = DummyModel_TPU()

params = SpectralFitting.all_parameter_symbols(model)
@test params == (:K, :a, :b)
free_symbols = SpectralFitting.free_parameter_symbols(model)
@test free_symbols == (:K, :a)
frozen_symbols = SpectralFitting.frozen_parameter_symbols(model)
@test frozen_symbols == (:b,)

# model + model

# SpectralFitting.freeze_parameter(model, :K)
