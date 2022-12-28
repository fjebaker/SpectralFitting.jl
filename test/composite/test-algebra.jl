using Test
using SpectralFitting

include("../dummies.jl")

add_model = DummyAdditive()
mul_model = DummyMultiplicative()

# allowed operations
cm = add_model + add_model
@test modelkind(cm) == SpectralFitting.Additive()
cm = mul_model * add_model
@test modelkind(cm) == SpectralFitting.Additive()
cm = mul_model * mul_model
@test modelkind(cm) == SpectralFitting.Multiplicative()

# don't care so much about the error message just the error
@test_throws "" add_model * add_model
@test_throws "" add_model * mul_model
@test_throws "" mul_model + mul_model


# check if parsing works okay too
model = mul_model * add_model
info = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
# model order should be always right then left 
for (i, check) in zip(info, ([:K, :a, :b], [:a, :b]))
    @test i.symbols == check
end
for (i, check) in zip(info, ([:K, :a], [:a]))
    @test i.free == check
end
for (i, check) in zip(info, ([:b], [:b]))
    @test i.frozen == check
end
