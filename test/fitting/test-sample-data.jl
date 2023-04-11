using SpectralFitting
using Test

# might have to download the model data first
SpectralFitting.download_model_data(XS_Laor)

function load_data(spec, mission)
    data = SpectralDataset(mission, spec)
    mask_bad_channels!(data)
    regroup!(data, data.spectrum.grouping)
    # divide by energy
    SpectralFitting.normalize_counts!(data)
    mask_domain!(data, <(0.8))
    mask_domain!(data, >(10.0))
    data
end

# path to the data directory
data1 = load_data(joinpath(testdir, "xmm/pn_spec_grp.fits"), XmmNewtonEPIC())

model = PhotoelectricAbsorption() * XS_PowerLaw() + XS_Laor()
display(model)

# construct the model and data problem
prob = FittingProblem(model, data1)

result = fit(prob, LevenbergMarquadt())

# these have been checked and are the same as XSPEC
@test result.χ2 ≈ 484.55 atol = 0.1
@test result.u ≈ [8.9215e-05, 6.5673, 0.010446, 1.8422, 0.13760] rtol = 1e-3
