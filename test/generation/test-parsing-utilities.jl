using Test
using SpectralFitting

include("../dummies.jl")
include("../fuzz.jl")

model = DummyAdditive()

params = SpectralFitting.all_parameter_symbols(model)
@test params == (:K, :a, :b)
free_symbols = SpectralFitting.free_parameter_symbols(model)
@test free_symbols == (:K, :a)
frozen_symbols = SpectralFitting.frozen_parameter_symbols(model)
@test frozen_symbols == (:b,)

# model info
info = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
@test info.symbols == [:K, :a, :b]
@test info.free == [:K, :a]
@test info.frozen == [:b]

# model introspetion
T = SpectralFitting.FunctionGeneration.model_T(typeof(model))
@test T === typeof(FitParam(1.0))

T = SpectralFitting.numbertype(model)
@test T === typeof(1.0)

# fuzz with all models we have
for model in FUZZ_ALL_MODELS
    _ = SpectralFitting.all_parameter_symbols(model)
    _ = SpectralFitting.free_parameter_symbols(model)
    _ = SpectralFitting.frozen_parameter_symbols(model)
    _ = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
    TT = SpectralFitting.FunctionGeneration.model_T(typeof(model))
    @test TT === typeof(FitParam(1.0))
    TT = SpectralFitting.numbertype(model)
    @test TT === typeof(1.0)
end


# introspection for composite models
cm = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())

@test_throws "" SpectralFitting.all_parameter_symbols(cm)
@test_throws "" SpectralFitting.free_parameter_symbols(cm)
@test_throws "" SpectralFitting.frozen_parameter_symbols(cm)

T = SpectralFitting.FunctionGeneration.model_T(typeof(cm))
@test T == typeof(FitParam(1.0))
T = SpectralFitting.numbertype(cm)
@test T == typeof(1.0)

info = SpectralFitting.FunctionGeneration.getinfo(typeof(cm))
@test length(info) == 3

# test the lenses we've assembled
model = cm
@test eval(info[1].lens) == model.right.right
@test eval(info[2].lens) == model.right.left
@test eval(info[3].lens) == model.left

# test lenses work on very deeply nested models 
model = DummyMultiplicative(a=FitParam(1.0)) * DummyMultiplicative(a=FitParam(2.0)) * (
    DummyAdditive() + DummyMultiplicative(a=FitParam(3.0)) * (
        DummyAdditive() + DummyMultiplicative(a=FitParam(4.0)) * (
            DummyAdditive() + DummyMultiplicative(a=FitParam(5.0)) * (
                DummyAdditive() + DummyMultiplicative(a=FitParam(6.0)) * (
                    DummyAdditive() + DummyMultiplicative(a=FitParam(7.0)) * (
                        DummyAdditive()
                    )
                )
            )
        )
    )
)
info = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
@test eval(info[1].lens) == model.right.right.right.right.right.right.right.right.right.right.right
@test eval(info[8].lens) == model.right.right.right.right.left
@test length(info) == 13

# order of parameters when we query for them is important
model = DummyAdditiveWithManyFrozen()
res = SpectralFitting.FunctionGeneration.all_parameter_symbols(typeof(model))
@test res == (:K, :a, :b, :c, :d, :e, :f, :g, :h, :i, :j)

res = SpectralFitting.FunctionGeneration.free_parameter_symbols(typeof(model))
@test res == (:K, :h)

res = SpectralFitting.FunctionGeneration.frozen_parameter_symbols(typeof(model))
@test res == (:a, :b, :c, :d, :e, :f, :g, :i, :j)


# test closure capture is working well
model = DummyMultiplicative() * DummyAdditiveTableModel()
ga = SpectralFitting.FunctionGeneration.assemble_aggregate_info(typeof(model), Float64)
lens = SpectralFitting.FunctionGeneration.assemble_closures(ga, typeof(model))
@test eval(lens[1]) == model.right.table

# test works well with multiple closure models
model = DummyMultiplicative() * (DummyMultiplicativeTableModel() * DummyAdditiveTableModel())
ga = SpectralFitting.FunctionGeneration.assemble_aggregate_info(typeof(model), Float64)
lens = SpectralFitting.FunctionGeneration.assemble_closures(ga, typeof(model))
@test eval(lens[1]) == model.right.right.table
@test eval(lens[2]) == model.right.left.table