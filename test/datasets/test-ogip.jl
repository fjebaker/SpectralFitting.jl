using Test
using SpectralFitting

xmm_config = StandardOGIPConfig(rmf_matrix_index = 2, rmf_energy_index = 3)
nustar_config = StandardOGIPConfig(rmf_matrix_index = 3, rmf_energy_index = 2)

xmm_rmf_path = joinpath(testdir, "xmm/pn.rmf")
nustar_rmf_path = joinpath(testdir, "nustar/nu60001047002A01_sr.rmf")

# test we can read in the different response matrices correctly
rmf_xmm = OGIP.read_rmf(xmm_rmf_path, xmm_config)
rmf_nustar = OGIP.read_rmf(nustar_rmf_path, nustar_config)

# make sure the correct number of values have been written in
@test count(!=(0), rmf_xmm.matrix) == 928061
@test count(!=(0), rmf_nustar.matrix) == 8536159

# test reading the ancillary files
xmm_arf_path = joinpath(testdir, "xmm/pn.arf")
nustar_arf_path = joinpath(testdir, "nustar/nu60001047002A01_sr.arf")

arf_xmm = OGIP.read_ancillary_response(xmm_arf_path, xmm_config)
arf_nustar = OGIP.read_ancillary_response(nustar_arf_path, nustar_config)

@test arf_xmm.bins_low ≈ rmf_xmm.bins_low atol = 1e-7
@test arf_nustar.bins_low ≈ rmf_nustar.bins_low atol = 1e-3

# test reading spectra
xmm_spec_path = joinpath(testdir, "xmm/pn_spec_grp.fits")
nustar_spec_path = joinpath(testdir, "nustar/nu60001047002A01_sr_grp_simple.pha")

spec_xmm = OGIP.read_spectrum(xmm_spec_path, xmm_config)
spec_nustar = OGIP.read_spectrum(nustar_spec_path, nustar_config)

@test spec_xmm.telescope == "XMM"
@test spec_nustar.telescope == "NuSTAR"

@test spec_xmm.poisson_statistics
@test spec_nustar.poisson_statistics
