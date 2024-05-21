using Test
using SpectralFitting
using SpectralFitting: Reflection

include("../utils.jl")
include("../dummies.jl")
include("../fuzz.jl")

M = DummyAdditive() |> typeof

info = Reflection.get_info(M, :model)

@test length(info.closure_symbols) == 0
@test length(info.generated_closure_symbols) == 0
@test info.symbols == [:K, :a, :b]

params = Reflection.get_parameter_symbols(M)
@test params == (:K, :a, :b)

# fuzz with all models we have
for model in FUZZ_ALL_MODELS
    M = typeof(model)
    _ = Reflection.get_parameter_symbols(M)
    _ = Reflection.get_info(M, :model)
end

# introspection for composite models
cm = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())
CM = typeof(cm)

syms = Reflection.get_parameter_symbols(CM)
@test syms == [:K_1, :a_1, :b_1, :K_2, :a_2, :b_2, :a_3, :b_3]

info = Reflection.get_info(CM, :cm)

# test the assembled lenses
@test eval(info.models[1][2].lens) == getfield(getfield(cm, :right), :right)
@test eval(info.models[2][2].lens) == getfield(getfield(cm, :right), :left)
@test eval(info.models[3][2].lens) == getfield(cm, :left)

# test lenses work on very deeply nested models 
cm =
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

info = Reflection.get_info(typeof(cm), :cm)
@test length(info.models) == 13

# order of parameters when we query for them is important
model = DummyAdditiveWithManyFrozen()
res = Reflection.get_parameter_symbols(typeof(model))
@test res == (:K, :a, :b, :c, :d, :e, :f, :g, :h, :i, :j)

# test closure capture is working well
model = DummyMultiplicative() * DummyAdditiveTableModel()
info = Reflection.get_info(typeof(model), :model)
info.models[1][2].closure_symbols == [:table]

model =
    (DummyMultiplicative() * DummyMultiplicative() * DummyAdditive()) +
    DummyMultiplicative() * DummyAdditive()
info = Reflection.get_info(typeof(model), :model)
expr = "$(info.model_expression)"
@test expr == "(m3 * m2) * a2 + m1 * a1"

destructure = [s => (i.model, info.model_symbols[s]) for (s, i) in info.models]

@test destructure == [
    :a1 => (DummyAdditive{FitParam{Float64}}, [:K_1, :a_1, :b_1]),
    :m1 => (DummyMultiplicative{FitParam{Float64}}, [:a_2, :b_2]),
    :a2 => (DummyAdditive{FitParam{Float64}}, [:K_2, :a_3, :b_3]),
    :m2 => (DummyMultiplicative{FitParam{Float64}}, [:a_4, :b_4]),
    :m3 => (DummyMultiplicative{FitParam{Float64}}, [:a_5, :b_5]),
]

model = DummyAdditive()
destructure = eval(Reflection.assemble_parameter_named_tuple(typeof(model)))
@test destructure == (; K = model.K, a = model.a, b = model.b)

model = DummyAdditive() + DummyAdditive(K = FitParam(3.0))
destructure = eval(Reflection.assemble_parameter_named_tuple(typeof(model)))
@test length(destructure) == 6

# remaking models
M = DummyAdditive() |> typeof
new = Reflection.reassemble_model(M, Vector{Float64})
@test "$new" == "(DummyAdditive){Float64}(parameters[1], parameters[2], parameters[3])"

M = DummyAdditiveTableModel() |> typeof
new = Reflection.reassemble_model(M, Vector{Float64})
