# simple binned dataset with no mission, background, or responses
# simple_spectrum just contains counts

function SpectralDataset(bins_low, bins_high, simple_spectrum, simple_errors)
    println("Creating simple spectral dataset")
    
    println("Creating spectrum")
    spectrum = Spectrum(
        collect(1:length(simple_spectrum)), # channels
        fill(1, size(simple_spectrum)), # quality
        fill(1, size(simple_spectrum)), # grouping
        simple_spectrum, # values
        "counts", # unit_string
        1.0, # exposure time
        1.0, # background scale
        1.0, # area scale
        SpectralFitting.ErrorStatistics.Gaussian, # assume Gaussian errors
        simple_errors, # errors
        0.0, # systematic error
        "no telescope", # telescope
        "no instrument", # instrument
    )
    println("Spectrum values are")
    println(spectrum.values)

    println("Creating basic metadata")
    meta = BasicMetadata("no path", "no path", "no telescope", "no instrument")

    println("Creating dataset")
    SpectralDataset(
        SpectralUnits.u"counts",
        meta,
        bins_low,
        bins_high,
        spectrum, # <-- should be Spectrum type
        missing, # no background
        missing, # no response
        missing, # no ancillary response
        BitVector(fill(true, size(spectrum.values)))
     )
end
