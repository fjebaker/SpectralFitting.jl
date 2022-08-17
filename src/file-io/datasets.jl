function SpectralDataset(
    meta::AbstractMetadata,
    rm::ResponseMatrix,
    counts,
    countserror,
    channels;
    T::Type = Float64,
)
    energy_bins_low, energy_bins_high = augmented_energy_channels(channels, rm)

    counts_vec::Vector{T} = T.(counts)
    countserr_vec::Vector{T} = T.(countserror)
    channels_vec::Vector{Int} = Int.(channels)

    SpectralDataset(
        energy_bins_low,
        energy_bins_high,
        counts_vec,
        countserr_vec,
        channels_vec,
        BitVector([true for _ in counts_vec]),
        meta,
        rm,
    )
end

function Base.propertynames(::SpectralDataset)
    (fieldnames(SpectralDataset)..., :energy_bin_widths)
end

function Base.getproperty(data::SpectralDataset, s::Symbol)
    if s in (:channels, :energy_bins_low, :energy_bins_high, :counts, :countserror)
        @views getfield(data, s)[data.mask]
    elseif s == :energy_bin_widths
        get_energy_bin_widths(data)
    else
        getfield(data, s)
    end
end

function get_energy_bin_widths(data::SpectralDataset)
    data.energy_bins_high .- data.energy_bins_low
end

# interface

function mask_bad_channels!(data::SpectralDataset)
    bad_indices = data.meta.quality .!= 0
    data.mask[bad_indices] .= false
    data
end

function _group_dataset(data::SpectralDataset{T,M}, grouping) where {T,M}
    indices = grouping_to_indices(grouping)
    N = length(indices) - 1
    @show N

    energy_low = zeros(T, N)
    energy_high = zeros(T, N)
    new_counts = zeros(T, N)
    new_errs = zeros(T, N)
    new_mask = zeros(Bool, N)

    um_energy_bins_low = getfield(data, :energy_bins_low)
    um_energy_bins_high = getfield(data, :energy_bins_high)
    um_counts = getfield(data, :counts)

    grouping_indices_callback(indices) do (i, index1, index2)
        energy_low[i] = um_energy_bins_low[index1]
        energy_high[i] = um_energy_bins_high[index2]

        selection = @views um_counts[index1:index2]
        new_counts[i] = sum(selection)
        new_errs[i] = Statistics.std(selection)

        # if not all are masked out, set true
        new_mask[i] = !all(==(false), data.mask[index1:index2])
    end

    rm = group_response(data.response, grouping)
    meta = group_meta(data.meta, data.mask, indices)

    SpectralDataset(
        energy_low,
        energy_high,
        new_counts,
        # todo: better error estimate?
        new_errs,
        collect(1:N),
        BitVector(new_mask),
        meta,
        rm,
    )
end

function set_max_energy!(data::SpectralDataset, emax)
    inds = data.energy_bins_high .≤ emax
    trim_dataset!(data::SpectralDataset, inds)
    count(!, inds)
end

function set_min_energy!(data::SpectralDataset, emin)
    inds = data.energy_bins_low .≥ emin
    trim_dataset!(data::SpectralDataset, inds)
    count(!, inds)
end

get_energy_bins(data::SpectralDataset{T,M}) where {T,M} = get_energy_bins(data, T)

get_energy_limits(x) = extrema(get_energy_bins(x))

# printing

function Base.show(io::IO, dataset::SpectralDataset{T,M}) where {T,M}
    print(io, "SpectralDataset[$(missiontrait(M))N=$(nrow(dataset.dataframe))]")
end

function Base.show(io::IO, ::MIME"text/plain", dataset::SpectralDataset{T,M}) where {T,M}
    e_min = Printf.@sprintf "%g" minimum(dataset.energy_bins_low)
    e_max = Printf.@sprintf "%g" maximum(dataset.energy_bins_high)

    rmf_e_min = Printf.@sprintf "%g" minimum(dataset.response.energy_bins_low)
    rmf_e_max = Printf.@sprintf "%g" maximum(dataset.response.energy_bins_high)

    rate_min = Printf.@sprintf "%g" minimum(dataset.counts)
    rate_max = Printf.@sprintf "%g" maximum(dataset.counts)

    descr = """SpectralDataset with $(length(dataset.channels)) populated channels:
       Mission: $(typeof(missiontrait(M)))
        . E min  : $(rpad(e_min, 12)) E max : $(e_max)
        . Counts : $(rpad(rate_min, 12)) to      $(rate_max)
       Instrument Response
        . $(length(dataset.response.channels)) RMF Channels
        . E min  : $(rpad(rmf_e_min, 12)) E max : $(rmf_e_max)
    """
    print(io, descr)
end
