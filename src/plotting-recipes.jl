# plotting
import MakieCore: @recipe
import Makie

@recipe(SpectralScatter) do scene
    Makie.Attributes(color = :black, positive_only = false)
end

function Makie.plot!(ss::SpectralScatter{<:Tuple{<:AbstractCommonSpectralDataset}})
    data = ss[1][]

    energy = (maxenergybins(data) .+ minenergybins(data)) ./ 2
    energy_error = energybinwidths(data) ./ 2

    counts = countbins(data)
    counts_error_low = counterrors(data)
    counts_error_high = deepcopy(counts_error_low)

    if ss.positive_only[]
        mask = counts .> 0

        energy = energy[mask]
        energy_error = energy_error[mask]
        counts = counts[mask]
        counts_error_low = counts_error_low[mask]
        counts_error_high = counts_error_high[mask]

        # make sure lower error bars > 0
        low_diff = counts .- counts_error_low
        # these tolerances must be hardcoded unfortunately, due to Float32s
        low_mask = low_diff .< 1e-8
        # and 1e-7 is best we can do on the y-axis
        counts_error_low[low_mask] .+= low_diff[low_mask] .- 1e-6
    end

    Makie.errorbars!(
        ss,
        energy,
        counts,
        counts_error_low,
        counts_error_high;
        color = ss.color,
    )
    Makie.errorbars!(ss, energy, counts, energy_error; color = ss.color, direction = :x)

    ss
end

@recipe(SpectralLine) do scene
    Makie.Attributes(color = :orange, positive_only = false)
end

function Makie.plot!(
    sl::SpectralLine{<:Tuple{<:AbstractVector{<:Number},<:AbstractVector{<:Number}}},
)
    energy = sl[1][]
    energy = @views(energy[1:end-1])
    counts = sl[2][]

    if sl.positive_only[]
        mask = energy .> 0
        energy = energy[mask]
        counts = counts[mask]
    end

    Makie.stairs!(sl, energy, counts; step = :post, color = sl.color)
end

function tick_formatter(values)
    map(values) do v
        if 10 ≥ v ≥ 1e-2
            "$v"
        else
            sup = round(Int, log10(v))
            s = Makie.UnicodeFun.to_superscript(sup)
            "10$s"
        end
    end
end

function makie_spectral_plot_axis(fig_loc; limits = (nothing, nothing), xticks = nothing)
    if isnothing(xticks)
        xticks = ([2, 5, 10, 100, 1000])
    end
    Makie.Axis(
        fig_loc,
        xscale = log10,
        yscale = log10,
        limits = limits,
        ytickformat = tick_formatter,
        xticks = xticks,
        xlabel = "Energy (keV)",
        ylabel = "Counts (s⁻¹ keV⁻¹)",
    )
end

function spectralplot(df::AbstractCommonSpectralDataset)
    fig = Makie.Figure()
    spectralplot(fig[1, 1], df)
    fig
end

function spectralplot(fig_loc, df::AbstractCommonSpectralDataset)
    limits = (getElims(df), (1e-6, maximum(countbins(df)) + 1))
    ax = makie_spectral_plot_axis(fig_loc; limits = limits)
    spectralscatter!(ax, df; positive_only = true)
    ax
end

export spectralscatter,
    spectralscatter!, makie_spectral_plot_axis, spectralplot, spectralline!
