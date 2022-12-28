using Test
using SpectralFitting


SpectralFitting.download_all_model_data()

ALL_XSPEC_MODELS = [
    XS_PowerLaw,
    XS_BlackBody,
    XS_BremsStrahlung,
    XS_Laor,
    XS_DiskLine,
    XS_KerrDisk,
    XS_KyrLine,
    XS_Gaussian
]

energy = collect(range(0.1, 100.0, 100))

for M in ALL_XSPEC_MODELS
    m = M()
    # can we invoke okay
    f = invokemodel(energy, m)
    @test all(.!isnan.(f))
    @test all(.!isinf.(f))
end