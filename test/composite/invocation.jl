using Test
using SpectralFitting

include("../dummies.jl")

model = DummyAdditive()

SpectralFitting.FunctionGeneration.generated_model_call!(Vector{Float64}, Vector{Float64}, model, Vector{Float64}, Vector{Float64})