using Test
using SpectralFitting

include("../dummies.jl")

dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
model = PowerLaw()

sim = simulate(model, dummy_data; seed = 42)
