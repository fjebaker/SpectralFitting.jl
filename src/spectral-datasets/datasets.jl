export mask_bad_channels!, mask_energy!, normalize_counts!

# utility constructor
function SpectralDataset(mission::AbstractMission, path; kwargs...)
    paths = read_OGIP_paths_from_file(path)
    SpectralDataset(mission, paths.spectrum, paths.response, paths.ancillary; kwargs...)
end

function SpectralDataset(
    units,
    spec::OGIP_GroupedEvents,
    rmf::OGIP_RMF,
    arf,
    background,
    meta::AbstractMetadata,
)
    rm = ResponseMatrix(rmf)
    ebins_low, ebins_high = augmented_energy_channels(spec.channels, rm)
    SpectralDataset(
        ebins_low,
        ebins_high,
        spec.values,
        spec.stat_error,
        units,
        meta,
        spec.header.poisson_error,
        rm,
        arf,
        background,
        spec.channels,
        spec.grouping,
        spec.quality,
        BitVector([true for _ in spec.values]),
        spec.header.exposure_time,
    )
end

# observation_id
function observation_id(::SpectralDataset{T,MetaType}) where {T,MetaType}
    error("Not implemented for mission $(MetaType)")
end

function observation_object(::SpectralDataset{T,MetaType}) where {T,MetaType}
    error("Not implemented for mission $(MetaType)")
end

function Base.propertynames(::SpectralDataset)
    (fieldnames(SpectralDataset)..., :counts, :rate, :energy_bin_widths)
end

unmasked(data::SpectralDataset, s) = getfield(data, s)

unmasked_counts(data::SpectralDataset{T,M,P,SpectralUnits._counts}, s) where {T,M,P} =
    getfield(data, s)

unmasked_counts(data::SpectralDataset{T,M,P,SpectralUnits._rate}, s) where {T,M,P} =
    getfield(data, s) .* data.exposure_time

unmasked_rate(data::SpectralDataset{T,M,P,SpectralUnits._counts}, s) where {T,M,P} =
    getfield(data, s) ./ data.exposure_time

unmasked_rate(data::SpectralDataset{T,M,P,SpectralUnits._rate}, s) where {T,M,P} =
    getfield(data, s)

function Base.getproperty(data::SpectralDataset, s::Symbol)
    if s in (:bins_low, :bins_high, :_errors, :_data)
        getfield(data, s)[data.mask]
    elseif s == :counts
        unmasked_counts(data, :_data)[data.mask]
    elseif s == :rate
        unmasked_rate(data, :_data)[data.mask]
    elseif s == :countserror
        unmasked_counts(data, :_errors)[data.mask]
    elseif s == :rateerror
        unmasked_rate(data, :_errors)[data.mask]
    elseif s == :channels
        channels = @views getfield(data, s)[data.mask]
        # channels must start at 0
        channels .- (minimum(channels))
    elseif s == :energy_bin_widths
        get_energy_bin_widths(data)
    else
        getfield(data, s)
    end
end

function normalize_counts!(data::SpectralDataset)
    ΔE = getfield(data, :bins_high) .- getfield(data, :bins_low)
    unmasked_counts(data, :_data) ./= ΔE
    unmasked_counts(data, :_errors) ./= ΔE
end

function fold_ancillary(data::SpectralDataset{T,M,P,U,A}) where {T,M,P,U,A}
    if A === Nothing
        data.response.matrix
    else
        data.ancillary.spec_response' .* data.response.matrix
    end
end

function get_energy_bin_widths(data::SpectralDataset)
    data.bins_high .- data.bins_low
end

# interface

function mask_bad_channels!(data::SpectralDataset)
    bad_indices = data.quality .!= 0
    data.mask[bad_indices] .= false
    data
end

function mask_energy!(data::SpectralDataset, cond)
    inds = cond.(unmasked(data, :bins_low))
    data.mask[inds] .= false
    data
end


function regroup(data::SpectralDataset{T,M,P,U}, grouping) where {T,M,P,U}
    indices = grouping_to_indices(grouping)
    N = length(indices) - 1

    energy_low = zeros(T, N)
    energy_high = zeros(T, N)
    new_data = zeros(T, N)
    new_errs = zeros(T, N)
    new_mask = zeros(Bool, N)

    um_bins_low = unmasked(data, :bins_low)
    um_bins_high = unmasked(data, :bins_high)
    um_data = unmasked(data, :_data)
    # um_errors = unmasked(data, :_errors)

    grouping_indices_callback(indices) do (i, index1, index2)
        energy_low[i] = um_bins_low[index1]
        energy_high[i] = um_bins_high[index2]

        selection = @views um_data[index1:index2]
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
    e_min = Printf.@sprintf "%g" minimum(data.bins_low)
    e_max = Printf.@sprintf "%g" maximum(data.bins_high)
    nbins = length(data.bins_low)

    rmf_e_min = Printf.@sprintf "%g" minimum(data.response.bins_low)
    rmf_e_max = Printf.@sprintf "%g" maximum(data.response.bins_high)

    rate_min = Printf.@sprintf "%g" minimum(data.counts)
    rate_max = Printf.@sprintf "%g" maximum(data.counts)
    exposure_time = Printf.@sprintf "%g" data.exposure_time

    is_grouped = all(==(1), data.grouping) ? "yes" : "no"

    has_background = isnothing(data.background) ? "no" : "yes"

    has_ancillary = isnothing(data.ancillary) ? "no" : "yes"

    descr = """SpectralDataset with $(length(data.channels)) populated channels:
       Object            : $(observation_object(data))
       Observation ID    : $(observation_id(data))
       Mission: $(typeof(missiontrait(M)))
        . Exposure time  : $(exposure_time) s
        . Bins           : $(nbins)
        . E (min/max)    : ($(e_min), $(e_max)) keV
        . Data (min/max) : ($(rate_min), $(rate_max)) $(data.units)
        . Grouped        : $(is_grouped)
        . Background     : $(has_background)
       Instrument Response
        . $(length(data.response.channels)) RMF Channels
        . E (min/max)    : ($(rmf_e_min), $(rmf_e_max)) keV
        . Anc            : $(has_ancillary)
    """
    print(io, descr)
end
