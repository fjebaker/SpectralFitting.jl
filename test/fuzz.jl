
# wont fuzz with XSPEC models because downloading data and also outside of our control
FUZZ_ALL_MODELS = [
    DummyAdditive(),
    DummyMultiplicative(),
    DummyMultiplicativeTableModel(),
    PowerLaw(),
    BlackBody(),
    PhotoelectricAbsorption(),
]
