using SpectralFitting
using Test

# might have to download the model data first
SpectralFitting.download_model_data(XS_Laor)
SpectralFitting.download_model_data(PhotoelectricAbsorption)

function prepare_data!(data, low, high)
    drop_bad_channels!(data)
    regroup!(data)
    normalize!(data)
    mask_energies!(data, low, high)
    data
end

# path to the data directory
data1 = SpectralFitting.XmmData(
    SpectralFitting.XmmEPIC(),
    joinpath(testdir, "xmm/pn_spec_grp.fits"),
)
prepare_data!(data1, 0.8, 10.0)

model = PhotoelectricAbsorption() * XS_PowerLaw() + XS_Laor()

# construct the model and data problem
prob = FittingProblem(model, data1)

result = fit(prob, LevenbergMarquadt())
# these have been checked and are the same as XSPEC
@test result.χ2 ≈ 484.88 atol = 0.1
@test result.u ≈ [8.9215e-05, 6.5673, 0.010446, 1.8422, 0.13760] rtol = 1e-3

data2 = NuStarData(joinpath(testdir, "nustar/nu60001047002A01_sr_grp_simple.pha"))
prepare_data!(data2, 2.0, 14.0)

prob = FittingProblem(model, data2)
result = fit(prob, LevenbergMarquadt())

@test result.χ2 ≈ 191.96312275498227 atol = 0.1
@test result.u ≈ [
    0.0002059507411466161,
    6.499205973061382,
    0.01989636801877023,
    1.9972614257180774,
    0.33187746370113763,
] rtol = 1e-3

# test joint
prob = FittingProblem(model => data1, model => data2)
result = fit(prob, LevenbergMarquadt())

@test result[1].χ2 ≈ 484.88 atol = 0.1
@test result[2].χ2 ≈ 191.96 atol = 0.1

# with binding
bind!(prob, :K_1)
result = fit(prob, LevenbergMarquadt())

@test result[1].χ2 ≈ 485.35 atol = 0.1
@test result[2].χ2 ≈ 203.19 atol = 0.1

# todo: with background subtraction
data1_nobkg = deepcopy(data1)
subtract_background!(data1_nobkg)

prob = FittingProblem(model, data1_nobkg)
result = fit(prob, LevenbergMarquadt())
# these have been checked and are the same as XSPEC
@test result.χ2 ≈ 496.6315394006663 atol = 0.1
@test result.u ≈ [9.2113e-05, 6.5597, 0.010478, 1.8483, 0.13960] rtol = 1e-3
