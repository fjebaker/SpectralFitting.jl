using Test
using SpectralFitting

# this model has been problematic so it gets its own explicit test case
model = XS_PhotoelectricAbsorption() * (XS_PowerLaw() + XS_KerrDisk())

energy = collect(range(0.2, 8.0, 1000))
fluxes = zeros(Float64, (length(energy) - 1, flux_count(model)))
flux = view(fluxes, :, 1)

# check model invokes okay
invokemodel!(fluxes, energy, model)
@test all(.!isnan.(flux))
@test all(.!isinf.(flux))

# check free parameters can be modified externally
check_flux = deepcopy(flux)
params = convert.(Float64, freeparameters(model))
invokemodel!(fluxes, energy, model, params)

@test isapprox.(check_flux, flux) |> all
