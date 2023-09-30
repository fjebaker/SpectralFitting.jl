struct SpectralDataPaths
    spectrum::Union{Missing,String}
    background::Union{Missing,String}
    response::Union{Missing,String}
    ancillary::Union{Missing,String}
end

function SpectralDataPaths(; spectrum = "", background = "", response = "", ancillary = "")
    SpectralDataPaths(
        spectrum == "" ? missing : spectrum,
        background == "" ? missing : background,
        response == "" ? missing : response,
        ancillary == "" ? missing : ancillary,
    )
end

function SpectralDataPaths(spec_path)
    background_path, response_path, ancillary_path =
        OGIP.read_paths_from_spectrum(spec_path)
    SpectralDataPaths(
        spectrum = spec_path,
        background = background_path,
        response = response_path,
        ancillary = ancillary_path,
    )
end

mutable struct SpectralData{T} <: AbstractDataset
    spectrum::Spectrum{T}
    response::ResponseMatrix{T}
    # background is optional
    background::Union{Missing,Spectrum{T}}
    # ancillary response is optionally, may also have already been folded into response
    ancillary::Union{Missing,AncillaryResponse{T}}

    channel_to_energy::Vector{T} # energy translated from the response channels
    domain::Vector{T} # domain fitted in models

    data_mask::BitVector
    domain_mask::BitVector
end

# constructor

SpectralData(paths::SpectralDataPaths, config::OGIP.AbstractOGIPConfig; kwargs...) =
    _dataset_from_ogip(paths, config; kwargs...)

function SpectralData(
    spectrum::Spectrum,
    response::ResponseMatrix;
    background = missing,
    ancillary = missing,
)
    domain = _make_domain_vector(spectrum, response)
    energy = _make_energy_vector(spectrum, response)
    data_mask = BitVector(fill(true, size(spectrum.data)))
    domain_mask = BitVector(fill(true, size(domain)))
    SpectralData(
        spectrum,
        response,
        background,
        ancillary,
        energy,
        domain,
        data_mask,
        domain_mask,
    )
end

supports_contiguosly_binned(::Type{<:SpectralData}) = true

function make_objective(layout::AbstractLayout, dataset::SpectralData)
    if dataset.spectrum.unit_string != "count / s keV"
        @warn "Spectrum is currently still in $(dataset.spectrum.unit_string). Most models fit in rate (count / s keV). Use `normalize!(dataset)` to ensure the dataset is in a standard format."
    end
    make_objective(layout, dataset.spectrum)[dataset.data_mask]
end

function make_objective_variance(layout::AbstractLayout, dataset::SpectralData)
    if dataset.spectrum.unit_string != "count / s keV"
        @warn "Spectrum is currently still in $(dataset.spectrum.unit_string). Most models fit in rate (count / s keV). Use `normalize!(dataset)` to ensure the dataset is in a standard format."
    end
    make_objective_variance(layout, dataset.spectrum)[dataset.data_mask]
end

make_domain(::ContiguouslyBinned, dataset::SpectralData) =
    dataset.domain[dataset.domain_mask]

restrict_domain!(dataset::SpectralData, low, high) =
    restrict_domain!(dataset, i -> high > i > low)
mask_energies!(dataset::SpectralData, low, high) =
    mask_energies!(dataset, i -> high > i > low)

function mask_energies!(dataset::SpectralData, condition)
    J = @. !condition(dataset.channel_to_energy)
    dataset.data_mask[@views(J[1:end-1])] .= false
    dataset
end

function restrict_domain!(dataset::SpectralData, condition)
    mask_energies!(dataset, condition)
    # I = @. !condition(dataset.domain)
    # dataset.domain_mask[I] .= false
    dataset
end

function objective_transformer(::ContiguouslyBinned, dataset::SpectralData)
    R = fold_ancillary(dataset.response, dataset.ancillary)[dataset.data_mask, :]
    ΔE = bin_widths(dataset)
    function _transformer!!(energy, flux)
        flux = R * flux
        @. flux = flux / ΔE
    end
    function _transformer!!(output, energy, flux)
        mul!(output, R, flux)
        @. output = output / ΔE
    end
    _transformer!!
end

bin_widths(dataset::SpectralData) = diff(dataset.channel_to_energy)[dataset.data_mask]
has_background(dataset::SpectralData) = !ismissing(dataset.background)
has_ancillary(dataset::SpectralData) = !ismissing(dataset.ancillary)

function drop_bad_channels!(dataset::SpectralData)
    indices = findall(==(BAD_QUALITY), dataset.spectrum.quality)
    drop_channels!(dataset, indices)
end

function drop_channels!(dataset::SpectralData, indices)
    drop_channels!(dataset.spectrum, indices)
    drop_channels!(dataset.response, indices)
    if has_background(dataset)
        drop_channels!(dataset.background, indices)
    end
    deleteat!(dataset.data_mask, indices)
    deleteat!(dataset.channel_to_energy, indices)
    dataset
end

spectrum_energy(dataset::SpectralData) =
    @views dataset.channel_to_energy[1:end-1][dataset.data_mask]
make_label(dataset::SpectralData) = "data"

function regroup!(dataset::SpectralData, grouping; safety_copy = false)
    grp::typeof(grouping) = if safety_copy
        copy(grouping)
    else
        grouping
    end

    itt = GroupingIterator(grp)
    last = first(itt)
    for i in itt
        dataset.channel_to_energy[i[1]] = dataset.channel_to_energy[i[2]]
        last = i
    end
    dataset.channel_to_energy[last[1]+1] = dataset.channel_to_energy[last[3]]

    if has_background(dataset)
        regroup!(dataset.background, grp)
    end
    regroup!(dataset.response, grp)
    regroup!(dataset.spectrum, grp)

    resize!(dataset.data_mask, length(itt))
    resize!(dataset.channel_to_energy, length(itt) + 1)
    # set everything to unmasked
    dataset.data_mask .= 1
    dataset
end

regroup!(dataset::SpectralData) = regroup!(dataset, dataset.spectrum.grouping)

function normalize!(dataset::SpectralData)
    ΔE = bin_widths(dataset)
    normalize!(dataset.spectrum)
    if dataset.spectrum.unit_string != "count / s keV"
        @. dataset.spectrum.data /= ΔE
        @. dataset.spectrum.errors /= ΔE
        dataset.spectrum.unit_string = "count / s keV"
    end
    if has_background(dataset)
        normalize!(dataset.background)
        if dataset.background.unit_string != "count / s keV"
            @. dataset.background.data /= ΔE
            @. dataset.background.errors /= ΔE
            dataset.background.unit_string = "count / s keV"
        end
    end
    dataset
end

# internal methods

function _subtract_background(dataset::SpectralData)
    # should all already be rates
    # errors added in quadrature
    # TODO: this needs fixing to propagate errors properly
    data_variance = @views dataset.spectrum.errors[dataset.data_mask] .^ 2
    background_variance = @views dataset.background.errors[dataset.data_mask] .^ 2
    _subtract_background(data_variance, background_variance, dataset)
end

function _subtract_background(data_variance, background_variance, dataset)
    aD = dataset.spectrum.area_scale
    bD = dataset.spectrum.background_scale
    aB = dataset.background.area_scale
    bB = dataset.background.background_scale
    @. (data_variance / aD) - (bD / bB) * (background_variance / aB)
end

function _dataset_from_ogip(paths::SpectralDataPaths, config::OGIP.AbstractOGIPConfig)
    spec = OGIP.read_spectrum(paths.spectrum, config)
    back = if !ismissing(paths.background)
        OGIP.read_background(paths.background, config)
    else
        @warn "No background file specified."
        missing
    end
    resp = if !ismissing(paths.response)
        OGIP.read_rmf(paths.response, config)
    else
        @warn "No response file specified."
        missing
    end
    ancillary = if !ismissing(paths.ancillary)
        OGIP.read_ancillary_response(paths.ancillary, config)
    else
        @warn "No ancillary file specified."
        missing
    end

    # convert everything to rates
    if spec.unit_string == "counts"
        spec.unit_string = "count / s"
        @. spec.data /= spec.exposure_time
        if !ismissing(spec.errors)
            @. spec.errors /= spec.exposure_time
        end
        if !ismissing(back) && back.unit_string == "counts"
            @. back.data /= back.exposure_time
        end
    end
    SpectralData(spec, resp; background = back, ancillary = ancillary)
end

function _make_domain_vector(::Spectrum, resp::ResponseMatrix{T}) where {T}
    domain = zeros(T, length(resp.bins_low) + 1)
    # todo: check these are indeed contiguous
    domain[2:end] .= resp.bins_high
    domain[1] = resp.bins_low[1]
    domain
end

function _make_energy_vector(spec::Spectrum, resp::ResponseMatrix{T}) where {T}
    augmented_energy_channels(
        spec.channels,
        resp.channels,
        resp.channel_bins_low,
        resp.channel_bins_high,
    )
end

macro _forward_SpectralData_api(args)
    if args.head !== :.
        error("Bad syntax")
    end
    T, field = args.args
    quote
        SpectralFitting.supports_contiguosly_binned(t::Type{<:$(T)}) = true
        SpectralFitting.make_domain(layout::SpectralFitting.AbstractLayout, t::$(T)) =
            SpectralFitting.make_domain(layout, getproperty(t, $(field)))
        SpectralFitting.make_domain_variance(
            layout::SpectralFitting.AbstractLayout,
            t::$(T),
        ) = SpectralFitting.make_domain_variance(layout, getproperty(t, $(field)))
        SpectralFitting.make_objective(layout::SpectralFitting.AbstractLayout, t::$(T)) =
            SpectralFitting.make_objective(layout, getproperty(t, $(field)))
        SpectralFitting.make_objective_variance(
            layout::SpectralFitting.AbstractLayout,
            t::$(T),
        ) = SpectralFitting.make_objective_variance(layout, getproperty(t, $(field)))
        SpectralFitting.objective_transformer(
            layout::SpectralFitting.AbstractLayout,
            t::$(T),
        ) = SpectralFitting.objective_transformer(layout, getproperty(t, $(field)))
        SpectralFitting.regroup!(t::$(T), args...) =
            SpectralFitting.regroup!(getproperty(t, $(field)), args...)
        SpectralFitting.restrict_domain!(t::$(T), args...) =
            SpectralFitting.restrict_domain!(getproperty(t, $(field)), args...)
        SpectralFitting.mask_energies!(t::$(T), args...) =
            SpectralFitting.mask_energies!(getproperty(t, $(field)), args...)
        SpectralFitting.drop_channels!(t::$(T), args...) =
            SpectralFitting.drop_channels!(getproperty(t, $(field)), args...)
        SpectralFitting.drop_bad_channels!(t::$(T)) =
            SpectralFitting.drop_bad_channels!(getproperty(t, $(field)))
        SpectralFitting.normalize!(t::$(T)) =
            SpectralFitting.normalize!(getproperty(t, $(field)))
        SpectralFitting.spectrum_energy(t::$(T)) =
            SpectralFitting.spectrum_energy(getproperty(t, $(field)))
        SpectralFitting.bin_widths(t::$(T)) =
            SpectralFitting.bin_widths(getproperty(t, $(field)))
    end |> esc
end

# printing utilities

function Base.show(io::IO, ::SpectralData{T}) where {T}
    print(io, "SpectralData")
end

function _printinfo(io, data::SpectralData{T}) where {T}
    domain = @views data.domain[data.domain_mask]
    descr = """SpectralData:
      . Masked channels     : $(count(==(false), data.data_mask))
      . Channel-Energy (CE) : $(extrema(data.channel_to_energy))
      . CE length           : $(length(data.channel_to_energy))
      . Masked domain       : $(count(==(false), data.domain_mask))
      . Domain length       : $(length(domain))
      . Domain (min/max)    : $(extrema(domain))
    """
    print(io, descr)
    print(io, "- Primary Spectrum:\n ")
    _printinfo(io, data.spectrum)
    if has_background(data)
        print(io, "- Background:\n ")
        _printinfo(io, data.background)
    else
        print(io, "- Background: missing\n")
    end
    print(io, "- Response:\n ")
    _printinfo(io, data.response)
    if has_ancillary(data)
        print(io, "- Ancillary:\n ")
        _printinfo(io, data.ancillary)
    else
        print(io, "- Ancillary: missing\n")
    end
end

export SpectralData,
    restrict_domain!, mask_energies!, drop_bad_channels!, drop_channels!, normalize!
