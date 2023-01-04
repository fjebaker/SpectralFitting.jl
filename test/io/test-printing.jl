using Test
using SpectralFitting

include("../dummies.jl")

# printing single
model = DummyAdditive()
buffer = IOBuffer()
Base.show(buffer, MIME"text/plain"(), model)
string = String(take!(buffer))
expected = """DummyAdditive
   K  => (1 ± 0.1)
   a  => (1 ± 0.1)
   b  => (5 ± 0.5)
"""
@test string == expected

# printing composite
model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())
buffer = IOBuffer()
Base.show(buffer, MIME"text/plain"(), model)
string = String(take!(buffer))
expected = """CompositeModel with 3 component models:
\e[36m     m1 * (a2 + a1)\e[0m
Model key and parameters:
   a1 => DummyAdditive
     K_1 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     a_1 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_1 ->  5\e[36m                     FROZEN\e[0m
   a2 => DummyAdditive
     K_2 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     a_2 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_2 ->  5\e[36m                     FROZEN\e[0m
   m1 => DummyMultiplicative
     a_3 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_3 ->  1\e[36m                     FROZEN\e[0m
"""
@test string == expected
