using RecipesBase

@recipe function _plotting_func(dataset::InjectiveData; data_layout = OneToOne())
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
    _yerr = if error_statistic(dataset) == ErrorStatistics.Poisson
        sqrt.(rateerror)
    else
        rateerror
    end
    yerr --> _yerr
    xerr --> bin_widths(dataset) ./ 2
    markerstrokecolor --> :auto
    if all(>(0), rate)
        yticks --> ([0.01, 0.1, 1, 10, 100], [0.01, 0.1, 1, 10, 100])
        yscale --> :log10
    end
    if all(>(0), rate)
        xticks --> ([1e-1, 1, 2, 5, 10, 20, 50, 100], [1e-1, 1, 2, 5, 10, 20, 50, 100])
        xscale --> :log10
    end
    xlabel --> "Energy (keV)"
    ylabel --> objective_units(dataset)
    label --> make_label(dataset)
    minorgrid --> true
    x = spectrum_energy(dataset)

    I = @. !isinf(x) && !isinf(rate)
    @views (x[I], rate[I])
end

plotting_domain(dataset::AbstractDataset) = spectrum_energy(dataset)
plotting_domain(dataset::InjectiveData) = dataset.domain

@recipe function _plotting_func(dataset::AbstractDataset, result::FittingResult)
    label --> "fit"
    seriestype --> :stepmid
    y = _f_objective(result.config)(result.config.domain, result.u)
    x = plotting_domain(dataset)
    if length(y) != length(x)
        error(
            "Domain mismatch. Are you sure you're plotting the result with the right dataset?",
        )
    end
    x, y
end

@recipe function _plotting_func(dataset::AbstractDataset, result::FittingResultSlice)
    label --> "fit"
    seriestype --> :stepmid
    y = invoke_result(result, result.u)
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
       !(typeof(r.args[1]) <: AbstractDataset) ||
       !(typeof(r.args[2]) <: AbstractFittingResult)
        error(
            "Ratio plots first argument must be `AbstractDataset` and second argument of type `AbstractFittingResult`.",
        )
    end

    data = r.args[1]
    x = plotting_domain(data)
    result = r.args[2] isa FittingResult ? r.args[2][1] : r.args[2]
    y = invoke_result(result, result.u)

    y_ratio = @. result.objective / y

    ylabel --> "Ratio [data / model]"
    xlabel --> "Energy (keV)"
    minorgrid --> true

    if (label == :auto)
        label = make_label(data)
    end

    @series begin
        linestyle --> :dash
        seriestype --> :hline
        label --> false
        color --> modelcolor
        [1.0]
    end

    @series begin
        markerstrokecolor --> datacolor
        label --> label
        seriestype --> :scatter
        markershape --> :none
        markersize --> 0.5
        yerror --> sqrt.(result.variance) ./ y
        xerror --> bin_widths(data) ./ 2
        x, y_ratio
    end
end
