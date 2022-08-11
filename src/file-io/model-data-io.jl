const MODEL_DATA_PATH = joinpath(LibXSPEC_jll.artifact_dir, "spectral",  "modelData")
MODEL_TO_MODEL_DATA_MAP = Dict{Symbol,Vector{String}}()
MODEL_DATA_PRESENT_CACHE = Dict{Symbol,Bool}()

function download_minimal_model_data()

end

function _register_model_data(M::Type{<:AbstractSpectralModel}, filenames::String...)
    s = Base.typename(M).name
    if s in keys(MODEL_TO_MODEL_DATA_MAP)
        push!(MODEL_TO_MODEL_DATA_MAP[s], filenames...)
    else
        MODEL_TO_MODEL_DATA_MAP[s] = collect(filenames)
    end
end

_is_model_data_downloaded(M::Type{<:AbstractSpectralModel}) = _is_model_data_downloaded(implementation(M), M)
_is_model_data_downloaded(::AbstractSpectralModelImplementation, _)::Bool = true
function _is_model_data_downloaded(::XSPECImplementation, M)::Bool
    s = Base.typename(M).name
    if get(MODEL_DATA_PRESENT_CACHE, s, false)
        return true
    else
        return _check_model_files_and_add_to_cache(s)
    end
end

function _check_model_files_and_add_to_cache(s::Symbol)
    # get the filenames we need
    filenames = get(MODEL_TO_MODEL_DATA_MAP, s, nothing)
    if !isnothing(filenames)
        if !all(i -> ispath(joinpath(MODEL_DATA_PATH, i)), filenames)
            return false
        end
    end
    # cache true
    MODEL_DATA_PRESENT_CACHE[s] = true
    return true
end

function _check_model_directory_present()
    if !ispath(MODEL_DATA_PATH)
        @warn "No model data directory found. Use `SpectralFitting.download_minimal_model_data()` to populate."
    end
end

function _download_from_archive(src, dest; progress=true, io::IO=Core.stdout)
    url = "http://www.star.bris.ac.uk/fbaker/XSPEC-model-data/$src"
    pg = if progress
        bar = MiniProgressBar(header="Downloading: $src", color=Base.info_color())
        start_progress(io, bar)
        (total, now) -> begin
            bar.max = total
            bar.current = now
            show_progress(io, bar)
        end
    else
        (_, _) -> nothing
    end
    Downloads.download(
        url, dest;
        progress = pg
    )
    end_progress(io, bar)
end

download_model_data(::M; kwargs...) where {M<:AbstractSpectralModel} = download_model_data(M; kwargs...)
function download_model_data(M::Type{<:AbstractSpectralModel}; kwargs...)
    if _is_model_data_downloaded(M)
        @warn "Model data for $(model_base_name(M)) is already downloaded."
        return false
    end

    # check model data directory exists
    if !ispath(MODEL_DATA_PATH)
        mkdir(MODEL_DATA_PATH)
    end

    s = Base.typename(M).name
    for src in MODEL_TO_MODEL_DATA_MAP[s]
        dest = joinpath(MODEL_DATA_PATH, src)
        if !ispath(dest)
            _download_from_archive(src, dest; kwargs...)
        end
    end
    true
end

export download_model_data
