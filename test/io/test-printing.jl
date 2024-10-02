using Test
using SpectralFitting

include("../utils.jl")
include("../dummies.jl")

function showstring(item)
    buffer = IOBuffer()
    Base.show(buffer, MIME"text/plain"(), item)
    String(take!(buffer))
end

# printing single
model = DummyAdditive()
string = showstring(model)
expected = """
┌ DummyAdditive
│      K ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      a ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      b ->  5  \e[36m                     FROZEN\e[0m
└ """
@test string == expected

# printing composite
model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())
# test destructuring works
# test the output string is correct
string = showstring(model)
expected = """┌ CompositeModel with 3 model components:
│ \e[36m     m1 * (a2 + a1)\e[0m
│ Model key and parameters:
│ \e[36m   a1\e[0m => \e[36mDummyAdditive\e[0m
│      K_1 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      a_1 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      b_1 ->  5  \e[36m                     FROZEN\e[0m
│ \e[36m   a2\e[0m => \e[36mDummyAdditive\e[0m
│      K_2 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      a_2 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      b_2 ->  5  \e[36m                     FROZEN\e[0m
│ \e[36m   m1\e[0m => \e[36mDummyMultiplicative\e[0m
│      a_3 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      b_3 ->  5  \e[36m                     FROZEN\e[0m
└ """
@test string == expected


# check string for more complex model
model =
    (DummyMultiplicative() * DummyMultiplicative() * DummyAdditive()) +
    DummyMultiplicative() * DummyAdditive()
string = showstring(model)
expected = """┌ CompositeModel with 5 model components:
│ \e[36m     (m3 * m2) * a2 + m1 * a1\e[0m
│ Model key and parameters:
│ \e[36m   a1\e[0m => \e[36mDummyAdditive\e[0m
│      K_1 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      a_1 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      b_1 ->  5  \e[36m                     FROZEN\e[0m
│ \e[36m   m1\e[0m => \e[36mDummyMultiplicative\e[0m
│      a_2 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      b_2 ->  5  \e[36m                     FROZEN\e[0m
│ \e[36m   a2\e[0m => \e[36mDummyAdditive\e[0m
│      K_2 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      a_3 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      b_3 ->  5  \e[36m                     FROZEN\e[0m
│ \e[36m   m2\e[0m => \e[36mDummyMultiplicative\e[0m
│      a_4 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      b_4 ->  5  \e[36m                     FROZEN\e[0m
│ \e[36m   m3\e[0m => \e[36mDummyMultiplicative\e[0m
│      a_5 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
│      b_5 ->  5  \e[36m                     FROZEN\e[0m
└ """
@test string == expected

# skip these for now
if false
    # test problem printing
    model = DummyMultiplicative() * DummyAdditive()
    data = make_dummy_dataset(x -> x)
    prob = FittingProblem(model, data)
    string = showstring(prob)
    expected = """┌ FittingProblem:
    │   Models:
    │     . CompositeModel[\e[36mDummyMultiplicative * DummyAdditive\e[0m]
    │   Data:
    │     . SpectralDataset[NoMission(),obs_id=noID]
    └ """
    @test string == expected

    # test results
    result = FittingResult([0.0, 1.0, 2.0, 3.0], 13.0, model, ones(Float64, 10), x -> x)
    string = showstring(result)
    expected = """┌ FittingResult:
    │     Model: CompositeModel[\e[36mDummyMultiplicative * DummyAdditive\e[0m]
    │     . u     : [0.0000, 1.0000, 2.0000, 3.0000]
    │     . χ²    : 13.000 
    └ """
    @test string == expected

    # mutli fitting results
    mr = MultiFittingResult((result, result, result))
    string = showstring(mr)
    expected = """┌ MultiFittingResult:
    │  FittingResult:
    │      Model: CompositeModel[\e[36mDummyMultiplicative * DummyAdditive\e[0m]
    │      . u     : [0.0000, 1.0000, 2.0000, 3.0000]
    │      . χ²    : 13.000 
    │  
    │  FittingResult:
    │      Model: CompositeModel[\e[36mDummyMultiplicative * DummyAdditive\e[0m]
    │      . u     : [0.0000, 1.0000, 2.0000, 3.0000]
    │      . χ²    : 13.000 
    │  
    │  FittingResult:
    │      Model: CompositeModel[\e[36mDummyMultiplicative * DummyAdditive\e[0m]
    │      . u     : [0.0000, 1.0000, 2.0000, 3.0000]
    │      . χ²    : 13.000 
    │  
    └ Σχ² = 39.000"""
    @test string == expected
end
