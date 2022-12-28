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
Symbols:    Models:              Params:
      a1 => DummyAdditive        (K_1, a_1, b_1) 
      a2 => DummyAdditive        (K_2, a_2, b_2) 
      m1 => DummyMultiplicative  (a_3, b_3) 
Symbols:     Values:
     K_1 =>  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     a_1 =>  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_1 =>  5\e[36m                     FROZEN\e[0m
     K_2 =>  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     a_2 =>  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_2 =>  5\e[36m                     FROZEN\e[0m
     a_3 =>  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_3 =>  5\e[36m                     FROZEN\e[0m
"""
@test string == expected
