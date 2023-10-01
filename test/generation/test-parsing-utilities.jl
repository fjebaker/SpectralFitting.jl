using Test
using SpectralFitting

include("../utils.jl")
include("../dummies.jl")
include("../fuzz.jl")

model = DummyAdditive()

params = SpectralFitting.all_parameter_symbols(model)
@test params == (:K, :a, :b)

# model info
info = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
@test info.symbols == [:K, :a, :b]

# model introspetion
T = SpectralFitting.FunctionGeneration.model_T(typeof(model))
@test T === typeof(FitParam(1.0))

T = SpectralFitting.numbertype(model)
@test T === typeof(1.0)

# fuzz with all models we have
for model in FUZZ_ALL_MODELS
    _ = SpectralFitting.all_parameter_symbols(model)
    _ = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
    TT = SpectralFitting.FunctionGeneration.model_T(typeof(model))
    @test TT === typeof(FitParam(1.0))
    TT = SpectralFitting.numbertype(model)
    @test TT === typeof(1.0)
end


# introspection for composite models
cm = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())
@test SpectralFitting.all_parameter_symbols(cm) ==
      (:K_1, :a_1, :b_1, :K_2, :a_2, :b_2, :a_3, :b_3)

# check the composite version
T = SpectralFitting.FunctionGeneration.model_T(typeof(cm))
@test T == typeof(FitParam(1.0))
T = SpectralFitting.numbertype(cm)
@test T == typeof(1.0)

info = SpectralFitting.FunctionGeneration.getinfo(typeof(cm))
@test length(info) == 3

# test the lenses we've assembled
model = cm
@test eval(info[1].lens) == getfield(getfield(model, :right), :right)
@test eval(info[2].lens) == getfield(getfield(model, :right), :left)
@test eval(info[3].lens) == getfield(model, :left)

# test lenses work on very deeply nested models 
model =
    DummyMultiplicative(a = FitParam(1.0)) *
    DummyMultiplicative(a = FitParam(2.0)) *
    (
        DummyAdditive() +
        DummyMultiplicative(a = FitParam(3.0)) * (
            DummyAdditive() +
            DummyMultiplicative(a = FitParam(4.0)) * (
                DummyAdditive() +
                DummyMultiplicative(a = FitParam(5.0)) * (
                    DummyAdditive() +
                    DummyMultiplicative(a = FitParam(6.0)) * (
                        DummyAdditive() +
                        DummyMultiplicative(a = FitParam(7.0)) * (DummyAdditive())
                    )
                )
            )
        )
    )
info = SpectralFitting.FunctionGeneration.getinfo(typeof(model))
@test length(info) == 13

# order of parameters when we query for them is important
model = DummyAdditiveWithManyFrozen()
res = SpectralFitting.FunctionGeneration.all_parameter_symbols(typeof(model))
@test res == (:K, :a, :b, :c, :d, :e, :f, :g, :h, :i, :j)

# test closure capture is working well
model = DummyMultiplicative() * DummyAdditiveTableModel()
ga = SpectralFitting.FunctionGeneration.assemble_aggregate_info(typeof(model), Float64)
lens = SpectralFitting.FunctionGeneration.assemble_closures(ga, typeof(model))
@test eval(lens[1]) == getfield(model, :right).table

# test works well with multiple closure models
model =
    DummyMultiplicative() * (DummyMultiplicativeTableModel() * DummyAdditiveTableModel())
ga = SpectralFitting.FunctionGeneration.assemble_aggregate_info(typeof(model), Float64)
lens = SpectralFitting.FunctionGeneration.assemble_closures(ga, typeof(model))

# test the flattening of composite models works
model =
    (DummyMultiplicative() * DummyMultiplicative() * DummyAdditive()) +
    DummyMultiplicative() * DummyAdditive()
t = SpectralFitting.FunctionGeneration._destructure_for_printing(typeof(model))
expr, info = eval(t)
@test expr == "(m3 * m2) * a2 + m1 * a1"
@test info == ((
    a1 = (DummyAdditive(), ("K_1", "a_1", "b_1")),
    m1 = (DummyMultiplicative(), ("a_2", "b_2")),
    a2 = (DummyAdditive(), ("K_2", "a_3", "b_3")),
    m2 = (DummyMultiplicative(), ("a_4", "b_4")),
    m3 = (DummyMultiplicative(), ("a_5", "b_5")),
))
