using Test
using SpectralFitting

testdir = get(ENV, "SF_TEST_SUITE_DATA", "./dev/spectral-fitting-test-suite/sample-data")

xmm_rmf_path = joinpath(testdir, "xmm/pn.rmf")
nustar_rmf_path = joinpath(testdir, "nustar/nu60001047002A01_sr.rmf")

# test we can read in the different dataset files correctly
data_xmm = SpectralFitting.OGIP.read_rmf(
    xmm_rmf_path,
    StandardOGIPConfig(rmf_matrix_index = 2, rmf_energy_index = 3),
)
data_xmm.matrix

length.(data_xmm.f_chan) |> maximum

data_nustar = SpectralFitting.OGIP.read_rmf(
    nustar_rmf_path,
    StandardOGIPConfig(rmf_matrix_index = 3, rmf_energy_index = 2),
)


rm = SpectralFitting.OGIP_RMF(xmm_rmf_path)
typeof(rm)
