using Test
using SpectralFitting

# test specific setup
struct DummyModel_TA{T,F} <: AbstractSpectralModel
    K::FitParam{T}
    a::FitParam{T}
    b::FitParam{T}
    function DummyModel_TA(;K = FitParam(1.0), a=FitParam(1.0), b=FitParam(5.0))
        new{SpectralFitting.parameter_type(K),SpectralFitting.FreeParameters{(:K,:a)}}(K, a, b)
    end
end
SpectralFitting.modelkind(::Type{<:DummyModel_TA}) = Additive()


model = DummyModel_TA()

info = SpectralFitting.assemble_aggregate_info(typeof(model))
@test length(info.parameters) == 3
@test length(info.closure_params) == 0
@test info.models == Type[typeof(model)]
@test info.maximum_flux_count == 1
@test length(info.statements) == 1

# composite additive models
new_model = model + model

info = SpectralFitting.assemble_aggregate_info(typeof(new_model))
@test length(info.parameters) == 6
@test length(info.closure_params) == 0
@test info.models == Type[typeof(model), typeof(model)]
@test info.maximum_flux_count == 2
@test length(info.statements) == 3
