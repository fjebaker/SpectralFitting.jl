using RecipesBase

@recipe function _plotting_func(d::SpectralDataset)
    seriestype := :scatter
    markersize := 0
    color := :black
    (rate, rateerror) = (d.rate ./ d.energy_bin_widths, d.rateerror ./ d.energy_bin_widths)
    yerr :=rateerror 
    xerr := d.energy_bin_widths ./ 2
    xscale := :log10
    xticks := ([1e-1, 1, 2, 5, 10, 20, 50, 100], [1e-1, 1, 2, 5, 10, 20, 50, 100])
    xlabel := "Energy (keV)"
    energy = (d.energy_bins_low .+ d.energy_bins_high) ./ 2
    (energy, rate)
end
