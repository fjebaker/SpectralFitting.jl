using Test
using SpectralFitting

include("../utils.jl")
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
# test destructuring works
out = @inferred SpectralFitting._destructure_for_printing(model)
@test out == (
    "m1 * (a2 + a1)",
    (
        a1 = (DummyAdditive(), ("K_1", "a_1", "b_1"), (true, true, false)),
        a2 = (DummyAdditive(), ("K_2", "a_2", "b_2"), (true, true, false)),
        m1 = (DummyMultiplicative(), ("a_3", "b_3"), (true, false)),
    ),
)
# test the output string is correct
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
     b_3 ->  5\e[36m                     FROZEN\e[0m
"""
@test string == expected


# check string for more complex model
model =
    (DummyMultiplicative() * DummyMultiplicative() * DummyAdditive()) +
    DummyMultiplicative() * DummyAdditive()
buffer = IOBuffer()
Base.show(buffer, MIME"text/plain"(), model)
string = String(take!(buffer))
expected = """CompositeModel with 5 component models:
\e[36m     (m3 * m2) * a2 + m1 * a1\e[0m
Model key and parameters:
   a1 => DummyAdditive
     K_1 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     a_1 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_1 ->  5\e[36m                     FROZEN\e[0m
   m1 => DummyMultiplicative
     a_2 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_2 ->  5\e[36m                     FROZEN\e[0m
   a2 => DummyAdditive
     K_2 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     a_3 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_3 ->  5\e[36m                     FROZEN\e[0m
   m2 => DummyMultiplicative
     a_4 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_4 ->  5\e[36m                     FROZEN\e[0m
   m3 => DummyMultiplicative
     a_5 ->  1 ± 0.1  ∈ [ 0, Inf ]\e[32m   FREE\e[0m
     b_5 ->  5\e[36m                     FROZEN\e[0m
"""
@test string == expected
