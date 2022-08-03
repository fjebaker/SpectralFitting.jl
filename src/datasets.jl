
# for future use: mission specific parsing
abstract type AbstractMissionTrait end
struct XMM_Newton <: AbstractMissionTrait end

abstract type AbstractSpectralDataset end

rmf(sd::AbstractSpectralDataset) = error("Not implemented for $(typeof(sd)).")
dataframe(sd::AbstractSpectralDataset) = error("Not implemented for $(typeof(sd)).")
countbins(sd::AbstractSpectralDataset) = error("Not implemented for $(typeof(sd)).")
counterrors(sd::AbstractSpectralDataset) = error("Not implemented for $(typeof(sd)).")
rmfenergybins(sd::AbstractSpectralDataset) = error("Not implemented for $(typeof(sd)).")
minenergybins(sd::AbstractSpectralDataset) = error("Not implemented for $(typeof(sd)).")
maxenergybins(sd::AbstractSpectralDataset) = error("Not implemented for $(typeof(sd)).")
energybinwidths(sd::AbstractSpectralDataset) = maxenergybins(sd) .- minenergybins(sd)
channels(sd::AbstractSpectralDataset) = error("Not implemented for $(typeof(sd)).")

dropbadchannels!(sd::AbstractSpectralDataset; kwargs...) =
    error("Not implemented for $(typeof(sd)).")
setEmax!(sd::AbstractSpectralDataset, emax) = error("Not implemented for $(typeof(sd)).")
setEmin!(sd::AbstractSpectralDataset, emin) = error("Not implemented for $(typeof(sd)).")
normalizecounts!(sd::AbstractSpectralDataset) = error("Not implemented for $(typeof(sd)).")
getElims(sd::AbstractSpectralDataset) =
    (minimum(minenergybins(sd)), maximum(maxenergybins(sd)))

# some aliases
response(sd::AbstractSpectralDataset) = rmf(sd)

abstract type AbstractCommonSpectralDataset <: AbstractSpectralDataset end

rmf(sd::AbstractCommonSpectralDataset) = sd.rmf
dataframe(sd::AbstractCommonSpectralDataset) = sd.dataframe
countbins(sd::AbstractCommonSpectralDataset) = sd.dataframe.RATE
counterrors(sd::AbstractCommonSpectralDataset) = sd.dataframe.STAT_ERR
rmfenergybins(sd::AbstractCommonSpectralDataset) = energybins(sd.rmf)
minenergybins(sd::AbstractCommonSpectralDataset) = sd.dataframe.E_MIN
maxenergybins(sd::AbstractCommonSpectralDataset) = sd.dataframe.E_MAX
energybinwidths(sd::AbstractCommonSpectralDataset) = sd.dataframe.E_DELTA
channels(sd::AbstractCommonSpectralDataset) = sd.dataframe.CHANNEL

dropbadchannels!(sd::AbstractCommonSpectralDataset; quality = 0) =
    sd.dataframe = sd.dataframe[sd.dataframe.QUALITY.==quality, :]
setEmax!(sd::AbstractCommonSpectralDataset, emax) =
    sd.dataframe = sd.dataframe[sd.dataframe.E_MAX.<emax, :]
setEmin!(sd::AbstractCommonSpectralDataset, emin) =
    sd.dataframe = sd.dataframe[sd.dataframe.E_MIN.>emin, :]

function normalizecounts!(sd::AbstractCommonSpectralDataset)
    if !sd.normalized
        sd.dataframe.RATE = sd.dataframe.RATE ./ sd.dataframe.E_DELTA
        sd.dataframe.STAT_ERR = sd.dataframe.STAT_ERR ./ sd.dataframe.E_DELTA
        sd.normalized = true
    end
    sd
end

function Base.show(io::IO, dataset::AbstractCommonSpectralDataset)
    print(io, "SpectralDataset[N=$(nrow(dataset.dataframe))]")
end

function Base.show(io::IO, ::MIME"text/plain", dataset::AbstractCommonSpectralDataset)
    e_min = Printf.@sprintf "%g" minimum(dataset.dataframe.E_MIN)
    e_max = Printf.@sprintf "%g" maximum(dataset.dataframe.E_MAX)

    rmf_e_min = Printf.@sprintf "%g" minimum(dataset.rmf.ebins.E_MIN)
    rmf_e_max = Printf.@sprintf "%g" maximum(dataset.rmf.ebins.E_MAX)

    rate_min = Printf.@sprintf "%g" minimum(dataset.dataframe.RATE)
    rate_max = Printf.@sprintf "%g" maximum(dataset.dataframe.RATE)

    descr = """SpectralDataset with $(nrow(dataset.dataframe)) populated channels:
        . E min : $(rpad(e_min, 12)) E max : $(e_max)
        . Rate  : $(rpad(rate_min, 12)) to      $(rate_max)
        . Rate is normalized ? $(dataset.normalized)
       Instrument Response
        . $(nrow(dataset.rmf.ebins)) RMF Channels
        . E min : $(rpad(rmf_e_min, 12)) E max : $(rmf_e_max)
    """
    print(io, descr)
end

mutable struct SpectralDataset{P,R} <: AbstractCommonSpectralDataset
    dataframe::P
    rmf::R
    normalized::Bool
    SpectralDataset(dataframe::P, rmf::R; normalized = false) where {P,R} =
        new{P,R}(dataframe, rmf, normalized)
end

export SpectralDataset,
    AbstractCommonSpectralDataset,
    AbstractSpectralDataset,
    rmf,
    data,
    countbins,
    counterrors,
    channels,
    dropbadchannels!,
    setEmin!,
    setEmax!,
    getElims,
    normalizecounts!,
    rmfenergybins,
    minenergybins,
    maxenergybins
