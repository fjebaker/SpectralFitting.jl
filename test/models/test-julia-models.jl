using Test
using SpectralFitting


ALL_JULIA_MODELS = [PowerLaw, BlackBody, PhotoelectricAbsorption]

energy = collect(range(0.1, 100.0, 100))

# test that all of the models actually are invocable and finite
for M in ALL_JULIA_MODELS
    m = M()
    outflux = invokemodel(energy, m)
    @test all(.!isnan.(outflux))
    @test all(.!isinf.(outflux))
end

