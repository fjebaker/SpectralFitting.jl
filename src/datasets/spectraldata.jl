struct SpectralDataPaths
    spectrum::Union{Nothing,String}
    background::Union{Nothing,String}
    response::Union{Nothing,String}
    ancillary::Union{Nothing,String}
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(paths::SpectralDataPaths))
    descr = """SpectralFilePaths:
      . Spectrum       : $(paths.spectrum)
      . Response       : $(paths.response)
      . Background     : $(paths.background)
      . Ancillary      : $(paths.ancillary)
    """
    print(io, descr)
end

function SpectralDataPaths(;
    spectrum = nothing,
    background = nothing,
    response = nothing,
    ancillary = nothing,
)
    SpectralDataPaths(spectrum, background, response, ancillary)
end

function SpectralDataPaths(
    spec_path;
    background = nothing,
    response = nothing,
    ancillary = nothing,
)
    background_path, response_path, ancillary_path =
        OGIP.read_paths_from_spectrum(spec_path)
    SpectralDataPaths(
        spectrum = spec_path,
        background = isnothing(background) ? background_path : background,
        response = isnothing(response) ? response_path : response,
        ancillary = isnothing(ancillary) ? ancillary_path : ancillary,
    )
end

mutable struct SpectralData{T} <: AbstractDataset
    spectrum::Spectrum{T}
    response::ResponseMatrix{T}
    # background is optional
    background::Union{Nothing,Spectrum{T}}
    # ancillary response is optionally, may also have already been folded into response
    ancillary::Union{Nothing,AncillaryResponse{T}}

    energy_low::Vector{T} # energy translated from the response channels
    energy_high::Vector{T} # energy translated from the response channels
    domain::Vector{T} # domain fitted in models

    data_mask::BitVector
end

# constructor

function SpectralData(paths::SpectralDataPaths; kwargs...)
    spec, resp, back, anc = _read_all_ogip(paths; kwargs...)

    # convert everything to rates
    if spec.units == u"counts"
        spec.units = u"counts / s"
        @. spec.data /= spec.exposure_time
        if !isnothing(spec.errors)
            @. spec.errors /= spec.exposure_time
        end
    end

    if !isnothing(back) && back.units == u"counts"
        back.units = u"counts / s"
        @. back.data /= back.exposure_time
        if !isnothing(back.errors)
            @. back.errors /= back.exposure_time
        end
    end

    SpectralData(spec, resp; background = back, ancillary = anc)
end

function SpectralData(
    spectrum::Spectrum,
    response::ResponseMatrix;
    # try to match the domains of the response matrix to the data
    match_domains = true,
    background = nothing,
    ancillary = nothing,
)
    domain = _make_domain_vector(spectrum, response)
    energy_low, energy_high = _make_energy_vector(spectrum, response)
    data_mask = BitVector(fill(true, size(spectrum.data)))
    data = SpectralData(
        spectrum,
        response,
        background,
        ancillary,
        energy_low,
        energy_high,
        domain,
        data_mask,
    )
    if !check_domains(data)
        if match_domains
            match_domains!(data)
        else
            @warn "The spectrum and response domains are unmatched. Use `match_domains!` to remedy. Results will assume you know what you're doing"
        end
    end
    return data
end

supports_contiguosly_binned(::Type{<:SpectralData}) = true

function check_units_warning(units)
    if units != u"counts / (s * keV)"
        @warn "Data is currently still in $(units). Most models fit in rate (count / (s * keV)). Use `normalize!(dataset)` to ensure the dataset is in a standard format."
    end
end

function make_objective(layout::AbstractDataLayout, dataset::SpectralData)
    check_units_warning(dataset.spectrum.units)
    make_objective(layout, dataset.spectrum)[dataset.data_mask]
end

function make_objective_variance(layout::AbstractDataLayout, dataset::SpectralData)
    check_units_warning(dataset.spectrum.units)
    make_objective_variance(layout, dataset.spectrum)[dataset.data_mask]
end

make_model_domain(::ContiguouslyBinned, dataset::SpectralData) = dataset.domain
make_output_domain(::ContiguouslyBinned, dataset::SpectralData) =
    folded_energy(dataset.response)

restrict_domain!(dataset::SpectralData, low, high) =
    restrict_domain!(dataset, i -> high > i > low)
mask_energies!(dataset::SpectralData, low, high) =
    mask_energies!(dataset, i -> high > i > low)

function mask_energies!(dataset::SpectralData, condition)
    J = @. !condition(dataset.energy_low) || !condition(dataset.energy_high)
    dataset.data_mask[J] .= false
    dataset
end

function restrict_domain!(dataset::SpectralData, condition)
    mask_energies!(dataset, condition)
    dataset
end

function _fold_transformer(T::Type, layout, R, ΔE, E)
    cache = DiffCache(construct_objective_cache(layout, T, length(E), 1))
    function _transformer!!(energy, flux)
        f = rebin_if_different_domains!(get_tmp(cache, flux), E, energy, flux)
        f = R * f
        @. f = f / ΔE
    end
    function _transformer!!(output, energy, flux)
        f = rebin_if_different_domains!(get_tmp(cache, flux), E, energy, flux)
        mul!(output, R, f)
        @. output = output / ΔE
    end
    _transformer!!
end

function objective_transformer(
    layout::ContiguouslyBinned,
    dataset::SpectralData{T},
) where {T}
    R_folded = if has_ancillary(dataset)
        sparse(fold_ancillary(dataset.response, dataset.ancillary))
    else
        dataset.response.matrix
    end
    R = R_folded[dataset.data_mask, :]
    ΔE = bin_widths(dataset)
    model_domain = response_energy(dataset.response)
    _fold_transformer(T, layout, R, ΔE, model_domain)
end


unmasked_bin_widths(dataset::SpectralData) = dataset.energy_high .- dataset.energy_low
bin_widths(dataset::SpectralData) = unmasked_bin_widths(dataset)[dataset.data_mask]
has_background(dataset::SpectralData) = !isnothing(dataset.background)
has_ancillary(dataset::SpectralData) = !isnothing(dataset.ancillary)

function drop_bad_channels!(dataset::SpectralData)
    indices = findall(!=(GOOD_QUALITY), dataset.spectrum.quality)
    drop_channels!(dataset, indices)
end

function drop_negative_channels!(dataset::SpectralData)
    indices = findall(<(0), dataset.spectrum.data)
    drop_channels!(dataset, indices)
end

function drop_channels!(dataset::SpectralData, indices)
    drop_channels!(dataset.spectrum, indices)
    drop_channels!(dataset.response, indices)
    if has_background(dataset)
        drop_channels!(dataset.background, indices)
    end
    deleteat!(dataset.data_mask, indices)
    deleteat!(dataset.energy_low, indices)
    deleteat!(dataset.energy_high, indices)
    length(indices)
end

spectrum_energy(dataset::SpectralData) =
    ((dataset.energy_low.+dataset.energy_high)./2)[dataset.data_mask]

function regroup!(dataset::SpectralData, grouping; safety_copy = false)
    grp::typeof(grouping) = if safety_copy
        copy(grouping)
    else
        grouping
    end

    itt = GroupingIterator(grp)
    for i in itt
        dataset.energy_low[i[1]] = dataset.energy_low[i[2]]
        dataset.energy_high[i[1]] = dataset.energy_high[i[3]]
    end

    if has_background(dataset)
        regroup!(dataset.background, grp)
    end
    regroup!(dataset.response, grp)
    regroup!(dataset.spectrum, grp)

    resize!(dataset.data_mask, length(itt))
    resize!(dataset.energy_low, length(itt))
    resize!(dataset.energy_high, length(itt))
    # set everything to unmasked
    dataset.data_mask .= 1
    dataset
end

regroup!(dataset::SpectralData) = regroup!(dataset, dataset.spectrum.grouping)

function normalize!(dataset::SpectralData)
    ΔE = bin_widths(dataset)
    normalize!(dataset.spectrum)
    if !(dataset.spectrum.units == u"counts / (s * keV)")
        @. dataset.spectrum.data /= ΔE
        @. dataset.spectrum.errors /= ΔE
        dataset.spectrum.units = u"counts / (s * keV)"
    end
    if has_background(dataset)
        normalize!(dataset.background)
        if !(dataset.background.units == u"counts / (s * keV)")
            @. dataset.background.data /= ΔE
            @. dataset.background.errors /= ΔE
            dataset.background.units = u"counts / (s * keV)"
        end
    end
    dataset
end

function subtract_background!(dataset::SpectralData)
    if !has_background(dataset)
        error("No background to subtract. Did you already subtract the background?")
    end
    subtract_background!(dataset.spectrum, dataset.background)
    dataset.background = nothing
    dataset
end

function set_domain!(dataset::SpectralData, domain)
    dataset.domain = domain
end

objective_units(data::SpectralData) = data.spectrum.units

error_statistic(data::SpectralData) = error_statistic(data.spectrum)

# internal methods

function rebin_if_different_domains!(output, data_domain, model_domain, input)
    if length(data_domain) == length(model_domain)
        @. output = input
    else
        interpolated_rebin!(output, data_domain, input, model_domain)
    end
    output
end

function _read_all_ogip(paths::SpectralDataPaths; forgiving = false)
    spec = if !isnothing(paths.spectrum)
        OGIP.read_spectrum(paths.spectrum)
    else
        if !forgiving
            @warn "No spectrum file found."
        end
        nothing
    end
    back = if !isnothing(paths.background)
        OGIP.read_background(paths.background)
    else
        if !forgiving
            @warn "No background file found."
        end
        nothing
    end
    resp = if !isnothing(paths.response)
        OGIP.read_rmf(paths.response)
    else
        if !forgiving
            throw(
                "No response file found in the header. Response must be specified with the keyword `response=PATH`.",
            )
        else
            nothing
        end
    end
    ancillary = if !isnothing(paths.ancillary)
        OGIP.read_ancillary_response(paths.ancillary)
    else
        if !forgiving
            @warn "No ancillary file found."
        end
        nothing
    end

    (spec, resp, back, ancillary)
end

function _make_domain_vector(::Spectrum, resp::ResponseMatrix{T}) where {T}
    domain = zeros(T, length(resp.bins_low) + 1)
    # todo: check these are indeed contiguous
    domain[2:end] .= resp.bins_high
    domain[1] = resp.bins_low[1]
    domain
end

function _make_energy_vector(spec::Spectrum, resp::ResponseMatrix{T}) where {T}
    full_domain = augmented_energy_channels(
        spec.channels,
        resp.channels,
        resp.channel_bins_low,
        resp.channel_bins_high,
    )
    high = full_domain[2:end]
    # full domain becomes the low
    resize!(full_domain, length(high))
    full_domain, high
end

function check_domains(data::SpectralData)
    (length(data.spectrum.channels) == length(data.response.channels)) &&
        (all(i -> i in data.response.channels, data.spectrum.channels))
end

function match_domains!(data::SpectralData)
    # drop parts of the response matrix that aren't in the spectrum
    I = filter(i -> i ∈ data.spectrum.channels, data.response.channels)
    data.response = ResponseMatrix(
        data.response.matrix[I, :],
        data.response.channels[I],
        data.response.channel_bins_low[I],
        data.response.channel_bins_high[I],
        data.response.bins_low,
        data.response.bins_high,
    )
end

macro _forward_SpectralData_api(args)
    if args.head !== :.
        error("Bad syntax")
    end
    T, field = args.args
    quote
        SpectralFitting.supports_contiguosly_binned(t::Type{<:$(T)}) = true
        SpectralFitting.make_output_domain(
            layout::SpectralFitting.AbstractDataLayout,
            t::$(T),
        ) = SpectralFitting.make_output_domain(layout, getfield(t, $(field)))
        SpectralFitting.make_model_domain(
            layout::SpectralFitting.AbstractDataLayout,
            t::$(T),
        ) = SpectralFitting.make_model_domain(layout, getfield(t, $(field)))
        SpectralFitting.make_domain_variance(
            layout::SpectralFitting.AbstractDataLayout,
            t::$(T),
        ) = SpectralFitting.make_domain_variance(layout, getfield(t, $(field)))
        SpectralFitting.make_objective(
            layout::SpectralFitting.AbstractDataLayout,
            t::$(T),
        ) = SpectralFitting.make_objective(layout, getfield(t, $(field)))
        SpectralFitting.make_objective_variance(
            layout::SpectralFitting.AbstractDataLayout,
            t::$(T),
        ) = SpectralFitting.make_objective_variance(layout, getfield(t, $(field)))
        SpectralFitting.objective_transformer(
            layout::SpectralFitting.AbstractDataLayout,
            t::$(T),
        ) = SpectralFitting.objective_transformer(layout, getfield(t, $(field)))
        SpectralFitting.regroup!(t::$(T), args...) =
            SpectralFitting.regroup!(getfield(t, $(field)), args...)
        SpectralFitting.restrict_domain!(t::$(T), args...) =
            SpectralFitting.restrict_domain!(getfield(t, $(field)), args...)
        SpectralFitting.mask_energies!(t::$(T), args...) =
            SpectralFitting.mask_energies!(getfield(t, $(field)), args...)
        SpectralFitting.drop_channels!(t::$(T), args...) =
            SpectralFitting.drop_channels!(getfield(t, $(field)), args...)
        SpectralFitting.drop_bad_channels!(t::$(T)) =
            SpectralFitting.drop_bad_channels!(getfield(t, $(field)))
        SpectralFitting.drop_negative_channels!(t::$(T)) =
            SpectralFitting.drop_negative_channels!(getfield(t, $(field)))
        SpectralFitting.normalize!(t::$(T)) =
            SpectralFitting.normalize!(getfield(t, $(field)))
        SpectralFitting.objective_units(t::$(T)) =
            SpectralFitting.objective_units(getfield(t, $(field)))
        SpectralFitting.spectrum_energy(t::$(T)) =
            SpectralFitting.spectrum_energy(getfield(t, $(field)))
        SpectralFitting.bin_widths(t::$(T)) =
            SpectralFitting.bin_widths(getfield(t, $(field)))
        SpectralFitting.subtract_background!(t::$(T), args...) =
            SpectralFitting.subtract_background!(getfield(t, $(field)), args...)
        SpectralFitting.set_domain!(t::$(T), args...) =
            SpectralFitting.set_domain!(getfield(t, $(field)), args...)
        SpectralFitting.error_statistic(t::$(T)) =
            SpectralFitting.error_statistic(getfield(t, $(field)))
    end |> esc
end


# printing utilities

function Base.show(io::IO, @nospecialize(data::SpectralData))
    print(io, "SpectralData[$(data.spectrum.telescope_name)]")
end

function _printinfo(io, data::SpectralData{T}) where {T}
    domain = @views data.domain
    ce_min = @views prettyfloat.(minimum(data.energy_low[data.data_mask]))
    ce_max = @views prettyfloat.(maximum(data.energy_high[data.data_mask]))
    dom_min, dom_max = @views prettyfloat.(extrema(domain))
    @views println(
        io,
        Crayons.Crayon(foreground = :cyan),
        "SpectralData",
        Crayons.Crayon(reset = true),
        " with ",
        Crayons.Crayon(foreground = :cyan),
        length(data.energy_low[data.data_mask]),
        Crayons.Crayon(reset = true),
        " active channels:",
    )
    descr = """  . Chn. E (min/max)    : ($ce_min, $ce_max)
      . Masked channels     : $(count(==(false), data.data_mask)) / $(length(data.data_mask))
      . Model domain size   : $(length(domain))
      . Domain (min/max)    : ($dom_min, $dom_max)
    """
    print(io, descr)

    print(
        io,
        Crayons.Crayon(foreground = :cyan),
        "Primary Spectrum:",
        Crayons.Crayon(reset = true),
        "\n ",
    )
    _printinfo(io, data.spectrum)

    print(
        io,
        Crayons.Crayon(foreground = :cyan),
        "Response:",
        Crayons.Crayon(reset = true),
        "\n ",
    )
    _printinfo(io, data.response)

    if has_background(data)
        print(
            io,
            Crayons.Crayon(foreground = :cyan),
            "Background: ",
            Crayons.Crayon(reset = true),
            "\n ",
        )
        _printinfo(io, data.background)
    else
        print(
            io,
            Crayons.Crayon(foreground = :dark_gray),
            "Background: nothing",
            Crayons.Crayon(reset = true),
            "\n",
        )
    end

    if has_ancillary(data)
        print(
            io,
            Crayons.Crayon(foreground = :cyan),
            "Ancillary:",
            Crayons.Crayon(reset = true),
            "\n ",
        )
        _printinfo(io, data.ancillary)
    else
        print(
            io,
            Crayons.Crayon(foreground = :dark_gray),
            "Ancillary: nothing",
            Crayons.Crayon(reset = true),
            "\n",
        )
    end
end

export SpectralData,
    restrict_domain!,
    mask_energies!,
    drop_bad_channels!,
    drop_negative_channels!,
    drop_channels!,
    normalize!,
    subtract_background!,
    set_domain!
