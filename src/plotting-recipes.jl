using RecipesBase

@recipe function _plotting_func(
    dataset::InjectiveData;
    data_layout = OneToOne(),
)
    seriestype --> :scatter
    markersize --> 1.0
    markershape --> :none
    markerstrokecolor --> :auto
    yerr -> dataset.codomain_variance
    xerr -> dataset.domain_variance
    label --> make_label(dataset)
    minorgrid --> true
    dataset.domain, dataset.codomain
end

@recipe function _plotting_func(
    dataset::AbstractDataset;
    data_layout = ContiguouslyBinned(),
)
    seriestype --> :scatter
    markersize --> 0.5
    markershape --> :none
    (rate, rateerror) = (
        make_objective(data_layout, dataset),
        make_objective_variance(data_layout, dataset),
    )
    yerr --> rateerror
    xerr --> bin_widths(dataset) ./ 2
    markerstrokecolor --> :auto
    xscale --> :log10
    yscale --> :log10
    yticks --> ([0.01, 0.1, 1, 10, 100], [0.01, 0.1, 1, 10, 100])
    xticks --> ([1e-1, 1, 2, 5, 10, 20, 50, 100], [1e-1, 1, 2, 5, 10, 20, 50, 100])
    xlabel --> "Energy (keV)"
    ylabel --> "counts s⁻¹ keV⁻¹"
    label --> make_label(dataset)
    minorgrid --> true
    x = spectrum_energy(dataset)
    (x, rate)
end

plotting_domain(dataset::AbstractDataset) = spectrum_energy(dataset)
plotting_domain(dataset::InjectiveData) = dataset.domain

@recipe function _plotting_func(dataset::AbstractDataset, result::FittingResult)
    label --> "fit"
    seriestype --> :stepmid
    y = _f_objective(result.config)(result.config.domain, result.u)
    x = plotting_domain(dataset)
    x, y
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
    energy = data.bins_low[data.mask]
    fieldnames(typeof(r))

    ylabel --> "Ratio [data / model]"
    xlabel --> "Energy (keV)"
    minorgrid --> true

    if (label == :auto)
        label = observation_id(data)
    end

    @series begin
        linestyle --> :dash
        seriestype --> :hline
        label --> false
        color --> modelcolor
        # energy, ones(length(energy))
        [1.0]
    end

    ratio_flux = target_vector(data) ./ model_flux
    @series begin
        markerstrokecolor --> datacolor
        label --> label
        seriestype --> :scatter
        markershape --> :none
        markersize --> 0.5
        yerror --> sqrt.(target_variance(data)) ./ model_flux
        xerror --> get_bin_widths(data) ./ 2
        energy, ratio_flux
    end
end
