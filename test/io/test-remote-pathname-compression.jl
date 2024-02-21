using Test
using SpectralFitting

info = SpectralFitting.ModelDataInfo("test.fits", "test.fits")
@test info.compression == SpectralFitting._NoCompression

info = SpectralFitting.ModelDataInfo("test.fits.gz", "test.fits")
@test info.compression == SpectralFitting._CompressedGzip

@test_throws AssertionError info =
    SpectralFitting.ModelDataInfo("test.fits.gz", "test.fits.gz")

@test SpectralFitting._trim_compression_filename("testfile.gz") == "testfile"
@test SpectralFitting._trim_compression_filename("testfile") == "testfile"
