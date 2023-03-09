using RecipesBase

@recipe function _plotting_func(d::SpectralDataset)
    seriestype := :scatter
    markersize := 0.5
    markershape := :none
    (rate, rateerror) = (d.rate, d.rateerror)
    yerr := rateerror
    xerr := d.energy_bin_widths ./ 2
    markerstrokecolor := :auto
    xscale := :log10
    yscale := :log10
    yticks := ([0.01, 0.1, 1, 10, 100], [0.01, 0.1, 1, 10, 100])
    xticks := ([1e-1, 1, 2, 5, 10, 20, 50, 100], [1e-1, 1, 2, 5, 10, 20, 50, 100])
    xlabel := "Energy (keV)"
    ylabel := "counts s⁻¹ keV⁻¹"
    label := observation_id(d)
    minorgrid := true
    energy = (d.bins_low .+ d.bins_high) ./ 2
    (energy, rate)
end

@recipe function _plotting_func(d::SimpleDataset)
    seriestype := :scatter
    xscale := :log10
    yscale := :log10
    markerstrokewidth := 0
    xlabel := "x ($(d.x_units))"
    ylabel := "y ($(d.y_units))"
    label := d.name
    if !isnothing(d.x_err)
        @views xerr := d.x_err[1:end-1]
    end
    if !isnothing(d.y_err)
        yerr := d.y_err
    end
    minorgrid := true
    @views (d.x[1:end-1], d.y)
end

@recipe function _plotting_func(r::FittingResult)
    label := "fit"
    r.x[1:end-1], r.folded_invoke(r.x, r.u)
end
@recipe function _plotting_func(x::AbstractVector, r::FittingResult)
    label := "fit"
    x, r.folded_invoke(r.x, r.u)
end

# ratio plots
@userplot RatioPlot
@recipe function _plotting_func(
    r::RatioPlot;
    datacolor = :auto,
    modelcolor = :auto,
    label = :auto,
)
    if length(r.args) != 2 ||
       !(typeof(r.args[1]) <: SpectralDataset) ||
       !(typeof(r.args[2]) <: AbstractVector || typeof(r.args[2]) <: FittingResult)
        error(
            "Ratio plots first argument must be `SpectralDataset` and second argument of type `AbstractVector`.",
        )
    end
    data = r.args[1]
    model_flux = if (typeof(r.args[2]) <: FittingResult)
        result = r.args[2]
        result.folded_invoke(result.x, result.u)
    else
        r.args[2]
    end
    energy = data.bins_low
    fieldnames(typeof(r))

    ylabel := "Ratio [data / model]"
    xlabel := "Energy (keV)"
    minorgrid := true

    if (label == :auto)
        label = observation_id(data)
    end

    @series begin
        linestyle := :dash
        seriestype := :hline
        label := false
        color := modelcolor
        # energy, ones(length(energy))
        [1.0]
    end

    ratio_flux = data.rate ./ model_flux
    @series begin
        markerstrokecolor := datacolor
        label := label
        seriestype := :scatter
        markershape := :none
        markersize := 0.5
        yerror := data.rateerror ./ model_flux
        xerror := data.energy_bin_widths ./ 2
        energy, ratio_flux
    end
end
