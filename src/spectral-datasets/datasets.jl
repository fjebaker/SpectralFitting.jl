export mask_bad_channels!, mask_energy!

function SpectralDataset(
    units,
    spec::OGIP_Spectrum,
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
        spec.poisson_error,
        rm,
        arf,
        background,
        spec.channels,
        spec.grouping,
        spec.quality,
        BitVector([true for _ in spec.values]),
        spec.exposure_time,
    )
end

function Base.propertynames(::SpectralDataset)
    (fieldnames(SpectralDataset)..., :counts, :rate, :energy_bin_widths)
end

unmasked_counts(data::SpectralDataset{T,M,P,SpectralUnits._counts}, s) where {T,M,P} =
    getfield(data, s)

unmasked_counts(data::SpectralDataset{T,M,P,SpectralUnits._rate}, s) where {T,M,P} =
    getfield(data, s) .* data.exposure_time

unmasked_rate(data::SpectralDataset{T,M,P,SpectralUnits._counts}, s) where {T,M,P} =
    getfield(data, s) ./ data.exposure_time

unmasked_rate(data::SpectralDataset{T,M,P,SpectralUnits._rate}, s) where {T,M,P} =
    getfield(data, s)

function Base.getproperty(data::SpectralDataset, s::Symbol)
    if s in (:energy_bins_low, :energy_bins_high, :_errors, :_data)
        @views getfield(data, s)[data.mask]
    elseif s == :counts
        @views unmasked_counts(data, :_data)[data.mask]
    elseif s == :rate
        @views unmasked_rate(data, :_data)[data.mask]
    elseif s == :countserror
        @views unmasked_counts(data, :_errors)[data.mask]
    elseif s == :rateerror
        @views unmasked_rate(data, :_errors)[data.mask]
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

function fold_ancillary(data::SpectralDataset{T,M,P,U,A}) where {T,M,P,U,A}
    if A === Nothing
        data.response.matrix
    else
        data.ancillary.spec_response' .* data.response.matrix
    end
end

function get_energy_bin_widths(data::SpectralDataset)
    data.energy_bins_high .- data.energy_bins_low
end

# interface

function mask_bad_channels!(data::SpectralDataset)
    bad_indices = data.quality .!= 0
    data.mask[bad_indices] .= false
    data
end

function mask_energy!(data::SpectralDataset, cond)
    inds = cond.(unmasked_low_energy_bins(data))
    data.mask[inds] .= false
    data
end

unmasked_low_energy_bins(data) = getfield(data, :energy_bins_low)
unmasked_high_energy_bins(data) = getfield(data, :energy_bins_high)

function regroup(data::SpectralDataset{T,M,P,U}, grouping) where {T,M,P,U}
    indices = grouping_to_indices(grouping)
    N = length(indices) - 1

    energy_low = zeros(T, N)
    energy_high = zeros(T, N)
    new_counts = zeros(T, N)
    new_errs = zeros(T, N)
    new_mask = zeros(Bool, N)

    um_energy_bins_low = unmasked_low_energy_bins(data)
    um_energy_bins_high = unmasked_high_energy_bins(data)
    um_counts = unmasked_counts(data)

    grouping_indices_callback(indices) do (i, index1, index2)
        energy_low[i] = um_energy_bins_low[index1]
        energy_high[i] = um_energy_bins_high[index2]

        selection = @views um_counts[index1:index2]
        new_counts[i] = sum(selection)
        new_errs[i] = sum(sqrt, selection)

        # if not all are masked out, set true
        new_mask[i] = !all(==(false), data.mask[index1:index2])
    end

    response = group_response_channels(data.response, grouping)
    quality = group_quality_vector(data.mask, data.quality, indices)

    SpectralDataset(
        data.units,
        energy_low,
        energy_high,
        new_counts,
        new_errs,
        data.meta,
        data.poisson_errors,
        response,
        data.ancillary,
        data.background,
        collect(1:N),
        [1 for _ in new_counts],
        quality,
        BitVector(new_mask),
        data.exposure_time,
    )
end

# printing

function Base.show(io::IO, dataset::SpectralDataset{T,M}) where {T,M}
    print(io, "SpectralDataset[$(missiontrait(M))N=$(nrow(dataset.dataframe))]")
end

function Base.show(io::IO, ::MIME"text/plain", data::SpectralDataset{T,M}) where {T,M}
    e_min = Printf.@sprintf "%g" minimum(data.energy_bins_low)
    e_max = Printf.@sprintf "%g" maximum(data.energy_bins_high)

    rmf_e_min = Printf.@sprintf "%g" minimum(data.response.energy_bins_low)
    rmf_e_max = Printf.@sprintf "%g" maximum(data.response.energy_bins_high)

    rate_min = Printf.@sprintf "%g" minimum(data.counts)
    rate_max = Printf.@sprintf "%g" maximum(data.counts)
    exposure_time = Printf.@sprintf "%g" data.exposure_time

    is_grouped = all(==(1), data.grouping) ? "yes" : "no"

    has_background = isnothing(data.background) ? "no" : "yes"

    has_ancillary = isnothing(data.ancillary) ? "no" : "yes"

    descr = """SpectralDataset with $(length(data.channels)) populated channels:
       Mission: $(typeof(missiontrait(M)))
        . Exposure time: $(exposure_time) s
        . E min      : $(rpad(e_min, 12)) E max : $(e_max)
        . Data       : $(rpad(rate_min, 12)) to      $(rate_max) $(data.units)
        . Grouped    : $(is_grouped)
        . Background : $(has_background)
       Instrument Response
        . $(length(data.response.channels)) RMF Channels
        . E min      : $(rpad(rmf_e_min, 12)) E max : $(rmf_e_max)
        . Anc        : $(has_ancillary)
    """
    print(io, descr)
end
