using Test
using SpectralFitting


SpectralFitting.download_all_model_data(verbose = false)

ALL_XSPEC_MODELS = [
    XS_PowerLaw,
    XS_BlackBody,
    XS_BremsStrahlung,
    XS_Laor,
    XS_DiskLine,
    XS_KerrDisk,
    XS_KyrLine,
    XS_Gaussian,
    XS_PhotoelectricAbsorption,
    XS_WarmAbsorption,
]

ALL_XSPEC_CONVOLUTIONAL = [XS_CalculateFlux]

energy = collect(range(0.1, 100.0, 100))

for M in ALL_XSPEC_MODELS
    m = M()
    # can we invoke okay
    f = invokemodel(energy, m)
    @test all(.!isnan.(f))
    @test all(.!isinf.(f))
end


# have to check convolutional models on a non-zero input flux
for M in ALL_XSPEC_CONVOLUTIONAL
    m = M()
    f = ones(Float64, length(energy) - 1)
    invokemodel!(f, energy, m)
    @test all(.!isnan.(f))
    @test all(.!isinf.(f))
end
