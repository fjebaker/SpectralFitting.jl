using RecipesBase

@recipe function _plotting_func(d::SpectralDataset)
    seriestype := :scatter
    markersize := 0
    color := :black
    (rate, rateerror) = (d.rate, d.rateerror)
    yerr := rateerror
    xerr := d.energy_bin_widths ./ 2
    xscale := :log10
    yticks := ([0.01, 0.1, 1, 10, 100], [0.01, 0.1, 1, 10, 100])
    xticks := ([1e-1, 1, 2, 5, 10, 20, 50, 100], [1e-1, 1, 2, 5, 10, 20, 50, 100])
    xlabel := "Energy (keV)"
    ylabel := "counts s⁻¹ keV⁻¹"
    minorgrid := true
    energy = (d.bins_low .+ d.bins_high) ./ 2
    (energy, rate)
end
