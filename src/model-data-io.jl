_model_data_storage_path = joinpath(LibXSPEC_jll.artifact_dir, "spectral", "modelData")
_model_to_data_map = Dict{Symbol,Vector{String}}()
_model_available_memoize_cache = Dict{Symbol,Bool}()

"""
    SpectralFitting.download_all_model_data()

Downloads all model data for the models currently registered with [`SpectralFitting.register_model_data`](@ref).
Calls [`SpectralFitting.download_model_data`](@ref) to perform the download.
"""
function download_all_model_data(; verbose = true)
    for s in keys(_model_to_data_map)
        download_model_data(s; verbose = verbose)
    end
end

"""
    SpectralFitting.register_model_data(M::Type{<:AbstractSpectralModel}, filenames::String...)
    SpectralFitting.register_model_data(s::Symbol, filenames::String...)

Register `filenames` as model data associated with the model given by type `M` or symbol `s`.
This function does not download any files, but rather adds the relevant filenames to a lookup which
[`SpectralFitting.download_model_data`](@ref) consults when invoked, and consequently model data is only downloaded when needed.

!!! note
    It is good practice to use this method immediately after defining a new model with [`@xspecmodel`](@ref)
    to register any required datafiles from the HEASoft source code, and therefore keep relevant information together.

# Example

```julia
# by type
register_model_data(XS_Laor, "ari.mod")
# by symbol
register_model_data(:XS_KyrLine, "KBHline01.fits")
```
"""
register_model_data(M::Type{<:AbstractSpectralModel}, filenames::String...) =
    register_model_data(Base.typename(M).name, filenames...)
function register_model_data(s::Symbol, filenames::String...)
    if s in keys(_model_to_data_map)
        push!(_model_to_data_map[s], filenames...)
    else
        _model_to_data_map[s] = collect(filenames)
    end
end

_is_model_data_downloaded(M::Type{<:AbstractSpectralModel}) =
    _is_model_data_downloaded(Base.typename(M).name)
function _is_model_data_downloaded(s::Symbol)::Bool
    if get(_model_available_memoize_cache, s, false)
        return true
    else
        return _check_model_files_and_add_to_cache(s)
    end
end

function _check_model_files_and_add_to_cache(s::Symbol)
    # get the filenames we need
    filenames = get(_model_to_data_map, s, nothing)
    if !isnothing(filenames)
        if !all(i -> ispath(joinpath(_model_data_storage_path, i)), filenames)
            return false
        end
    end
    # cache true
    _model_available_memoize_cache[s] = true
    return true
end

function _check_model_directory_present()
    if !ispath(_model_data_storage_path)
        @warn "No model data directory found. Use `SpectralFitting.download_all_model_data()` to populate."
    end
end

function _download_from_archive(
    src,
    dest;
    progress = true,
    io::IO = Core.stdout,
    model_source_url = "https://www.star.bris.ac.uk/fergus/XSPEC-model-data/",
)
    url = "$model_source_url/$src"
    pg = if progress
        bar = MiniProgressBar(header = "Downloading: $src", color = Base.info_color())
        start_progress(io, bar)
        (total, now) -> begin
            bar.max = total
            bar.current = now
            show_progress(io, bar)
        end
    else
        (_, _) -> nothing
    end
    Downloads.download(url, dest; progress = pg)
    end_progress(io, bar)
end

"""
    SpectralFitting.download_model_data(model::AbstractSpectralModel; kwargs...)
    SpectralFitting.download_model_data(M::Type{<:AbstractSpectralModel}; kwargs...)
    SpectralFitting.download_model_data(s::Symbol; kwargs...)

Downloads the model data for a model specified either by `model`, type `M`, or symbol `s`. Datafiles
associated with a specific model may be registered using [`SpectralFitting.register_model_data`](@ref).
The download is currently unconfigurable, but permits slight control via a number of keyword arguments:

- `progress::Bool = true`

Display a progress bar for the download.

- `model_source_url::String = "http://www.star.bris.ac.uk/fbaker/XSPEC-model-data"`

The source URL used to download the model data.

All standard XSPEC spectral model data is currently being hosted on the University of Bristol
astrophysics servers, and should be persistently available to anyone.
"""
download_model_data(M::Type{<:AbstractSpectralModel}; kwargs...) =
    download_model_data(Base.typename(M).name; kwargs...)
download_model_data(::M; kwargs...) where {M<:AbstractSpectralModel} =
    download_model_data(Base.typename(M).name; kwargs...)
function download_model_data(s::Symbol; verbose = true, kwargs...)
    _infolog(s) =
        if verbose
            @info s
        end
    if _is_model_data_downloaded(s)
        _infolog("Model data for $(s) is already downloaded.")
        return nothing
    end

    # check model data directory exists
    if !ispath(_model_data_storage_path)
        # else make it
        mkdir(_model_data_storage_path)
    end

    _infolog("Checking model data for $s:")
    for src in _model_to_data_map[s]
        dest = joinpath(_model_data_storage_path, src)
        if !ispath(dest)
            _download_from_archive(src, dest; kwargs...)
            _infolog("$src downloaded")
        end
    end
    _infolog("All requisite model data for $s downloaded.")
    return nothing
end

function ensure_model_data(M::Type)
    if !_is_model_data_downloaded(M)
        @warn "Model data for $(model_base_name(M)) is not present!\nRequisite model data may be fetched with `SpectralFitting.download_model_data($(model_base_name(M)))`."
        error("Missing data.")
    end
end

function load_and_unpack_model_data(M)
    files = _model_to_data_map[Base.typename(M).name]
    # load all of the data files
    contents = map(files) do f
        path = joinpath(_model_data_storage_path, f)
        load(path)
    end
    contents
end

function get_model_data(M::Type)
    ensure_model_data(M)
    load_and_unpack_model_data(M)
end
