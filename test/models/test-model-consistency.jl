using Test
using SpectralFitting


# ensure PowerLaw and XS_PowerLaw behave exactly the same
energy = collect(range(0.1, 100.0, 20))
for a in [2.0, 0.3, 1.33, 1.0]
    pl_f1 = invokemodel(energy, PowerLaw(K = 1.0, a = a))
    pl_f2 = invokemodel(energy, XS_PowerLaw(K = 1.0, a = a))
    @test isapprox.(pl_f1 .- pl_f2, 0.0, atol = 1e-4) |> all
end

@ciskip begin # models with downloaded data
    energy = collect(range(0.1, 100.0, 20))

    for ηH in [0.0, 0.1, 1.0, 2.0]
        # current interpolation insn't good above 2.0
        phabs_f1 = invokemodel(energy, PhotoelectricAbsorption(ηH = ηH))
        phabs_f2 = invokemodel(energy, XS_PhotoelectricAbsorption(ηH = ηH))
        @test isapprox.(phabs_f1 .- phabs_f2, 0.0, atol = 1e-1) |> all
    end
end
