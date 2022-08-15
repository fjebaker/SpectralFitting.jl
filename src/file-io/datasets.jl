export AbstractSpectralDatasetMeta,
    trim_dataset_meta!,
    SpectralDataset,
    trim_dataset!,
    set_max_energy!,
    set_min_energy!,
    normalize_counts!,
    get_energy_limits

# abstract supertype
abstract type AbstractSpectralDatasetMeta end

# interface
trim_dataset_meta!(::AbstractSpectralDatasetMeta, inds) = nothing

# concrete type
mutable struct SpectralDataset{M,R,T,C}
    meta::M
    response::R
    low_energy_bins::T
    high_energy_bins::T
    counts::T
    countserror::T
    channels::C
    normalized::Bool
end

# get_response(sd::SpectralDataset) = sd.response
# get_low_energy_bins(sd::SpectralDataset) = sd.low_energy_bins
# get_high_energy_bins(sd::SpectralDataset) = sd.high_energy_bins
# get_counts(sd::SpectralDataset) = sd.counts
# get_countserror(sd::SpectralDataset) = sd.countserror
# get_channels(sd::SpectralDataset) = sd.channels

# interface

function trim_dataset!(sd::SpectralDataset, inds)
    sd.low_energy_bins = sd.low_energy_bins[inds]
    sd.high_energy_bins = sd.high_energy_bins[inds]
    sd.counts = sd.counts[inds]
    sd.countserror = sd.countserror[inds]
    sd.channels = sd.channels[inds]
    trim_dataset_meta!(sd.meta, inds)
end
function set_max_energy!(sd::SpectralDataset, emax)
    inds = sd.high_energy_bins .≤ emax
    trim_dataset!(sd::SpectralDataset, inds)
    count(!, inds)
end
function set_min_energy!(sd::SpectralDataset, emin)
    inds = sd.low_energy_bins .≥ emin
    trim_dataset!(sd::SpectralDataset, inds)
    count(!, inds)
end

# additional
function normalize_counts!(sd::SpectralDataset)
    if !sd.normalized
        ΔE = sd.high_energy_bins .- sd.low_energy_bins
        sd.counts ./= ΔE
        sd.countserror ./= ΔE
        sd.normalized = true
    end
end

get_energy_bins(sd::SpectralDataset{M,R,T}) where {M,R,T} = get_energy_bins(sd, T)

function get_energy_bins(x, T::Type)
    energy = zeros(T, length(x.low_energy_bins) + 1)
    energy[1:end-1] .= x.low_energy_bins
    energy[end] = x.high_energy_bins[end]
    energy
end

get_energy_limits(x) = extrema(get_energy_bins(x))

# printing

function Base.show(io::IO, dataset::SpectralDataset{M}) where {M}
    print(io, "SpectralDataset[$(missiontrait(M))N=$(nrow(dataset.dataframe))]")
end

function Base.show(io::IO, ::MIME"text/plain", dataset::SpectralDataset{M}) where {M}
    e_min = Printf.@sprintf "%g" minimum(dataset.low_energy_bins)
    e_max = Printf.@sprintf "%g" maximum(dataset.high_energy_bins)

    rmf_e_min = Printf.@sprintf "%g" minimum(dataset.response.low_energy_bins)
    rmf_e_max = Printf.@sprintf "%g" maximum(dataset.response.high_energy_bins)

    rate_min = Printf.@sprintf "%g" minimum(dataset.counts)
    rate_max = Printf.@sprintf "%g" maximum(dataset.counts)

    descr = """SpectralDataset with $(length(dataset.channels)) populated channels:
       Mission: $(typeof(missiontrait(M)))
        . E min : $(rpad(e_min, 12)) E max : $(e_max)
        . Rate  : $(rpad(rate_min, 12)) to      $(rate_max)
        . Rate is normalized ? $(dataset.normalized)
       Instrument Response
        . $(length(dataset.response.channels)) RMF Channels
        . E min : $(rpad(rmf_e_min, 12)) E max : $(rmf_e_max)
    """
    print(io, descr)
end
