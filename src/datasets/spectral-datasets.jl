export mask_bad_channels!, normalize_counts!

# utility constructor
function SpectralDataset(mission::AbstractMission, path; kwargs...)
    paths = path_assembler(mission)(path)
    SpectralDataset(mission, paths; kwargs...)
end

function dataset_from_ogip(paths, config, meta)
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
        @. spec.values /= spec.exposure_time
        if !ismissing(spec.errors)
            @. spec.errors /= spec.exposure_time
        end
        if !ismissing(back) && back.unit_string == "counts"
            @. back.values /= back.exposure_time
        end
    end

    blow, bhigh = augmented_energy_channels(
        spec.channels,
        resp.channels,
        resp.channel_bins_low,
        resp.channel_bins_high,
    )
    units = SpectralUnits.u"counts / s"
    mask = BitVector(fill(true, size(spec.values)))
    SpectralDataset(units, meta, blow, bhigh, spec, back, resp, ancillary, mask)
end

# observation_id
function observation_id(::SpectralDataset{T,MetaType}) where {T,MetaType}
    error("Not implemented for mission $(MetaType)")
end

function observation_object(::SpectralDataset{T,MetaType}) where {T,MetaType}
    error("Not implemented for mission $(MetaType)")
end

# incase anyone every knew what the old api was like
# https://docs.julialang.org/en/v1/manual/style-guide/#Prefer-exported-methods-over-direct-field-access
get_units(::SpectralDataset{T,M,Units}) where {T,M,Units} = Units
function get_rate(s::SpectralDataset)
    rate =
        get_units(s) <: SpectralUnits._counts ?
        s.spectrum.values / s.spectrum.exposure_time : s.spectrum.values
    rate[s.mask]
end
function get_counts(s::SpectralDataset)
    counts =
        get_units(s) <: SpectralUnits._counts ? s.spectrum.values :
        s.spectrum.values * s.spectrum.exposure_time
    counts[s.mask]
end
function get_rate_variance(s::SpectralDataset)
    err =
        get_units(s) <: SpectralUnits._counts ?
        s.spectrum.errors / s.spectrum.exposure_time : s.spectrum.errors
    err[s.mask] .^ 2
end
function get_count_variance(s::SpectralDataset)
    err =
        get_units(s) <: SpectralUnits._counts ? s.spectrum.errors :
        s.spectrum.errors * s.spectrum.exposure_time
    err[s.mask] .^ 2
end
function get_channels(s::SpectralDataset)
    if has_background(s)
        s.response.channels
    else
        error("No channels defined")
    end
end
function get_bin_widths(s::SpectralDataset)
    low, high = get_bins(s)
    @. high - low
end
function get_bins(s::SpectralDataset)
    s.bins_low[s.mask], s.bins_high[s.mask]
end
function get_exposure_time(s::SpectralDataset)
    s.spectrum.exposure_time
end
function get_response(s::SpectralDataset{T}) where {T}
    if !has_response(s)
        # TODO: better errors
        error("No RMF defined")
    else
        s.response
    end
end
function get_ancillary(s::SpectralDataset{T}) where {T}
    if !has_ancillary(s)
        # TODO: better errors
        error("No ancillary response defined")
    else
        s.ancillary
    end
end
has_response(s::SpectralDataset) = !ismissing(s.response)
has_ancillary(s::SpectralDataset) = !ismissing(s.ancillary)
has_background(s::SpectralDataset) = !ismissing(s.background)

target_vector(data::SpectralDataset) = get_rate(data)
target_variance(data::SpectralDataset) = get_rate_variance(data)
domain_vector(data::SpectralDataset) = domain_vector(data.response)
function _lazy_folded_invokemodel(model::AbstractSpectralModel, data::SpectralDataset)
    ΔE = get_bin_widths(data)
    # pre-mask the response matrix to ensure channel out corresponds to the active data points
    R = fold_ancillary(data)[data.mask, :]
    # pre-allocate the output 
    wrapped = (energy, params) -> begin
        flux = invokemodel(energy, model, params)
        flux = (R * flux)
        @. flux = flux / ΔE
    end
    wrapped
end

function normalize_counts!(data::SpectralDataset)
    ΔE = get_bin_widths(data)
    @. data.y /= ΔE
    @. data.y_err /= ΔE
end

function fold_ancillary(data::SpectralDataset)
    if !has_response(data)
        return LinearAlgebra.I
    end
    response = get_response(data)
    ancillary = if has_ancillary(data)
        get_ancillary(data)
    else
        return response.matrix
    end
    ancillary.spec_response' .* response.matrix
end

function mask_bad_channels!(data::SpectralDataset)
    bad_indices = data.spectrum.quality .!= 0
    data.mask[bad_indices] .= false
    data
end

function mask_domain!(data::SpectralDataset, f)
    to_mask = @. f(data.bins_low) | f(data.bins_high)
    data.mask[to_mask] .= false
    data
end

function regroup(data::SpectralDataset{U,M,T}, grouping) where {U,M,T}
    indices = grouping_to_indices(grouping)
    N = length(indices) - 1

    energy_low = zeros(T, N)
    energy_high = zeros(T, N)
    new_data = zeros(T, N)
    new_errs = zeros(T, N)
    new_mask = zeros(Bool, N)

    old_bins_low, old_bins_high = get_bins(data)
    old_data = data.y
    # um_errors = unmasked(data, :_errors)

    grouping_indices_callback(indices) do (i, index1, index2)
        energy_low[i] = old_bins_low[index1]
        energy_high[i] = old_bins_high[index2]

        selection = @views old_data[index1:index2]
        new_data[i] = sum(selection)
        new_errs[i] = count_error(sum(selection), 1.0)

        # if not all are masked out, set true
        new_mask[i] = !all(==(false), data.mask[index1:index2])
    end

    response = group_response_channels(data.response, grouping)
    quality = group_quality_vector(data.mask, data.quality, indices)

    SpectralDataset(
        energy_low,
        energy_high,
        new_data,
        new_errs,
        data.units,
        data.meta,
        data.poisson_errors,
        response,
        data.ancillary,
        data.background,
        collect(1:N),
        [1 for _ in new_data],
        quality,
        BitVector(new_mask),
        data.exposure_time,
    )
end

# printing

function Base.show(io::IO, dataset::SpectralDataset{T,M}) where {T,M}
    print(io, "SpectralDataset[$(missiontrait(M)),obs_id=$(observation_id(dataset))]")
end

function Base.show(io::IO, ::MIME"text/plain", data::SpectralDataset{T,M}) where {T,M}
    buff = IOBuffer()
    _printinfo(buff, data)
    s = String(take!(buff))
    print(io, encapsulate(s))
end

function _printinfo(io, data::SpectralDataset{T,M}) where {T,M}
    blow, bhigh = get_bins(data)
    e_min = Printf.@sprintf "%g" minimum(blow)
    e_max = Printf.@sprintf "%g" maximum(bhigh)
    nbins = length(data.bins_low)

    rmf_e_min = Printf.@sprintf "%g" minimum(data.response.bins_low)
    rmf_e_max = Printf.@sprintf "%g" maximum(data.response.bins_high)

    y = target_vector(data)
    rate_min = Printf.@sprintf "%g" minimum(y)
    rate_max = Printf.@sprintf "%g" maximum(y)
    exposure_time = Printf.@sprintf "%g" data.spectrum.exposure_time

    is_grouped = all(==(1), data.spectrum.grouping) ? "yes" : "no"

    has_anc = has_ancillary(data) ? "yes" : "no"

    descr = """SpectralDataset with $(length(data.spectrum.channels)) populated channels:
       Object            : $(observation_object(data))
       Observation ID    : $(observation_id(data))
       Mission: $(typeof(missiontrait(M)))
        . Exposure time  : $(exposure_time) s
        . Bins           : $(nbins)
        . E (min/max)    : ($(e_min), $(e_max)) keV
        . Data (min/max) : ($(rate_min), $(rate_max)) $(data.units)
        . Grouped        : $(is_grouped)
       Instrument Response
        . $(length(data.response.channels)) RMF Channels
        . E (min/max)    : ($(rmf_e_min), $(rmf_e_max)) keV
        . Ancillary      : $(has_anc)
    """
    print(io, descr)
    desc_background = if has_background(data)
        bg_min = Printf.@sprintf "%g" minimum(data.background.values)
        bg_max = Printf.@sprintf "%g" maximum(data.background.values)
        """   Background
            . Exposure time  : ($(bg_min), $(bg_max)) $(data.units)
            . Data (min/max) :
        """
    else
        """   No Background
        """
    end
    print(io, desc_background)
end
