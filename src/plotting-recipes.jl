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
    _yerr = sqrt.(rateerror)
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

# residual plots
# note: multiple datasets require repeated calls to this function (write a wrapper later)
@userplot ResidualPlot
@recipe function _plotting_fun(
    r::ResidualPlot,
    datacolor = :auto,
    modelcolor = :auto,
    residualcolor = :auto,
    label = :auto,
)
    # check that the function has been passed one dataset and one fit result
    if length(r.args) != 2 ||
       !(typeof(r.args[1]) <: AbstractDataset) ||
       !(typeof(r.args[2]) <: AbstractFittingResult)
        error(
            "Ratio plots first argument must be `AbstractDataset` and second argument of type `AbstractFittingResult`.",
        )
    end

    println("Debug: Creating a residual plot")

    data = r.args[1]
    x = plotting_domain(data)
    # at the moment I don't understand why the following line is necessary
    # I would assume result = r.args[2] which might be of type FittingResultSlice
    result = r.args[2] isa FittingResult ? r.args[2][1] : r.args[2]
    y = invoke_result(result, result.u)

    # residual is the difference between the model and the data in units of "sigma" so the error bars have size 1
    # this assumes we have statistics such that sigma = sqrt(variance) - should probably make this more statistically neutral
    yerr = sqrt.(result.variance)
    y_residual = @. (result.objective - y) / yerr
    # is this the best way to ensure y_residual_error has the same type as y_residual, or should it just be fixed at Float64?
    y_residual_error = ones(eltype(y_residual), length(y_residual))

    minorgrid --> true

    if (label == :auto)
        label = make_label(data)
    end

    # layout --> @layout [grid(2, 1, heights=[0.7 ,0.3]), margin=0mm]
    layout --> @layout [
        top{0.75h}
        bottom{0.25h}
    ]
    margins --> (0, :mm)

    # logarithmic x-axis (might want to let this be an option)
    xscale --> :log10
    filtered_array = filter(x -> x != 0, result.objective)
    min_data_value = 0.8 * minimum(filtered_array)
    max_data_value = 1.2 * maximum(filtered_array)
    filtered_array = filter(x -> x != 0, y)
    min_model_value = 0.8 * minimum(filtered_array)
    max_model_value = 1.2 * maximum(filtered_array)
    min_non_zero_value = min(min_data_value, min_model_value)
    max_non_zero_value = max(max_data_value, max_model_value)

    # plot the data
    @series begin
        subplot --> 1
        yscale --> :log10
        yrange --> [min_non_zero_value, max_non_zero_value]
        markerstrokecolor --> datacolor
        label --> label
        seriestype --> :scatter
        markershape --> :none
        markersize --> 0.5
        yerror --> yerr
        xerror --> bin_widths(data) ./ 2
        x, result.objective
    end

    # plot the model - need to fix this
    @series begin
        subplot --> 1
        yscale --> :log10
        yrange --> [min_non_zero_value, max_non_zero_value]
        xticks --> nothing
        ylabel --> "Flux (units)"
        markerstrokecolor --> modelcolor
        label --> :none
        seriestype --> :stepmid
        markershape --> :none
        markersize --> 0.5
        xerror --> bin_widths(data) ./ 2
        x, y
    end

    # plot the residuals
    @series begin
        subplot --> 2
        yscale --> :identity
        xticks --> true
        xlabel --> "Energy (keV)"
        ylabel --> "(Data - model)/error"
        markerstrokecolor --> residualcolor
        label --> :none
        seriestype --> :scatter
        markershape --> :none
        markersize --> 0.5
        yerror --> y_residual_error
        xerror --> bin_widths(data) ./ 2
        x, y_residual
    end

    # zero line
    @series begin
        subplot --> 2
        yscale --> :identity
        xticks --> true
        xlabel --> "Energy (keV)"
        ylabel --> "(Data - model)/error"
        linestyle --> :dash
        seriestype --> :hline
        label --> false
        # not currently specifying a colour
        [0.0]
    end
end
