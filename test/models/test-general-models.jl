using Test
using SpectralFitting

# this file tests for specific model constructions
# and anything that produces errors should be added to this
# test set

include("../dummies.jl")
include("../utils.jl")

@ciskip begin
    model = PhotoelectricAbsorption() * (PowerLaw() + BlackBody())

    energy = collect(range(0.2, 8.0, 1000))
    fluxes = zeros(Float64, (length(energy) - 1, objective_cache_count(model)))
    flux = view(fluxes, :, 1)

    # check model invokes okay
    invokemodel!(fluxes, energy, model)
    @test all(.!isnan.(flux))
    @test all(.!isinf.(flux))

    # check free parameters can be modified externally
    check_flux = deepcopy(flux)
    params = convert.(Float64, freeparameters(model))
    invokemodel!(fluxes, energy, model, params)
    # there should be no difference
    @test isapprox.(check_flux, flux) |> all
end

model = DummyMultiplicative() * DummyMultiplicative() * (DummyAdditive() + DummyAdditive())

checker_all(model)
