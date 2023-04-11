using Test
using SpectralFitting

xmm_path = joinpath(testdir, "xmm/pn_spec_grp.fits")

nustar_path = joinpath(testdir, "nustar/nu60001047002A01_sr_grp_simple.pha")

SpectralDataset(XmmNewtonEPIC(), xmm_path)
