using SpectralFitting
using Test

# might have to download the model data first
SpectralFitting.download_model_data(XS_Laor)
SpectralFitting.download_model_data(PhotoelectricAbsorption)

function _prepare_data(data)
    drop_bad_channels!(data)
    regroup!(data)
    normalize!(data)
    mask_energies!(data, 0.8, 10.0)
    data
end

# path to the data directory
data1 =
    SpectralFitting.XMMData(
        SpectralFitting.XmmEPIC(),
        joinpath(testdir, "xmm/pn_spec_grp.fits"),
    ) |> _prepare_data

model = PhotoelectricAbsorption() * XS_PowerLaw() + XS_Laor()

# construct the model and data problem
prob = FittingProblem(model, data1)

result = fit(prob, LevenbergMarquadt())
# these have been checked and are the same as XSPEC
@test result.χ2 ≈ 484.88 atol = 0.1
@test result.u ≈ [8.9215e-05, 6.5673, 0.010446, 1.8422, 0.13760] rtol = 1e-3

# todo: with background subtraction
# subtract_background!(data1)
# prob = FittingProblem(model, data1)
# result = fit(prob, LevenbergMarquadt())
# # these have been checked and are the same as XSPEC
# @test result.χ2 ≈ 496.09 atol = 0.1
# @test result.u ≈ [9.181e-5, 6.55961, 0.0104753, 1.848017, 0.139510] rtol = 1e-3
