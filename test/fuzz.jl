
# wont fuzz with XSPEC models because downloading data and also outside of our control
FUZZ_ALL_MODELS = [
    DummyAdditive(),
    DummyMultiplicative(),
    DummyMultiplicativeTableModel(),
    PowerLaw(),
    BlackBody(),
]

# has data requirements, so skip on the CI
if get(ENV, "CI", false) == false
    push!(FUZZ_ALL_MODELS, PhotoelectricAbsorption())
end
