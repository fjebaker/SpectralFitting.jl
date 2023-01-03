include("utils.jl")

# wont fuzz with XSPEC models because downloading data and also outside of our control
FUZZ_ALL_MODELS = [
    DummyAdditive(),
    DummyMultiplicative(),
    DummyMultiplicativeTableModel(),
    PowerLaw(),
    BlackBody(),
]

# has data requirements, so skip on the CI
@ciskip push!(FUZZ_ALL_MODELS, PhotoelectricAbsorption())
