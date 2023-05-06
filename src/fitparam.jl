"""
Utility structure for tracking which parameters are free or frozen.
"""
struct FreeParameters{V} end

function fit_param_default_error(val)
    # 10 % error
    round(abs(0.1 * val), sigdigits = 1)
end

# concrete types
mutable struct FitParam{T}
    value::T
    error::T

    lower_limit::T
    upper_limit::T

    frozen::Bool

    FitParam(
        val::T;
        lower_limit = T(0.0),
        upper_limit = T(Inf),
        error = fit_param_default_error(val),
        frozen = false,
    ) where {T} = new{T}(val, error, lower_limit, upper_limit, frozen)
end

struct ParameterCache{T}
    fitparams::Vector{FitParam{T}}
    params::Vector{T}
end

function Base.iterate(cache::ParameterCache, i)
    if i > length(cache.fitparams)
        return nothing 
    end
    return cache.fitparams[i], i + 1
end

function Base.iterate(cache::ParameterCache)
    i::Int = 1
    Base.iterate(cache, i)
end

function update_free!(cache::ParameterCache, free)
    i::Int = 1
    N = length(free)
    for (j, fp) in enumerate(cache)
        if isfree(fp)
            if i > N
                error("Too few free parameters given")
            end
            cache.params[j] = free[i]
            i += 1
        end
    end
    if N >= i
        error("Too many free parameters given")
    end
    cache
end

function update_free_and_frozen!(cache::ParameterCache, free, frozen)
    if length(cache.params) != length(free) + length(frozen)
        error("Too many/few free/frozen parameters given")
    end
    i1::Int = 1
    i2::Int = 1
    for (j, fp) in enumerate(cache)
        if isfree(fp)
            cache.params[j] = free[i1]
            i1 += 1 
        else
            cache.params[j] = frozen[i2]
            i2 += 1
        end
    end
    cache
end

# interface
set_value!(f::FitParam{T}, val::T) where {T} = f.value = val
set_error!(f::FitParam{T}, val::T) where {T} = f.error = val
get_value(f::FitParam) = f.value
function set!(f::FitParam, o::FitParam)
    f.value = o.value
    f.lower_limit = o.lower_limit
    f.upper_limit = o.upper_limit
    f.error = o.error
    f
end
# edge case
get_value(x::Number) = x

isfrozen(f::FitParam) = f.frozen
isfree(f::FitParam) = !isfrozen(f)

get_error(f::FitParam) = f.error
get_upperlimit(f::FitParam) = f.upper_limit
get_lowerlimit(f::FitParam) = f.lower_limit

Base.isapprox(f1::FitParam, f2::FitParam; kwargs...) =
    isapprox(f1.value, f2.value; kwargs...)
Base.:(==)(f1::FitParam, f2::FitParam) = f1.value == f2.value
Base.convert(T::Type{<:Number}, f::FitParam) = convert(T, f.value)

parameter_type(::Type{FitParam{T}}) where {T} = T
parameter_type(::T) where {T<:FitParam} = parameter_type(T)

function get_info_tuple(f::FitParam)
    s1 = Printf.@sprintf "%.3g" get_value(f)
    s2 = Printf.@sprintf "%.3g" get_error(f)
    s3 = Printf.@sprintf "%.3g" get_lowerlimit(f)
    s4 = Printf.@sprintf "%.3g" get_upperlimit(f)
    (s1, s2, s3, s4)
end
# todo: edge case that should be avoided
# but currently needed when printing composite models
get_info_tuple(n::Number) = (Printf.@sprintf("%.3g", n), "0", "0", "0")

function print_info(io::IO, f::FitParam)
    v, e, lb, ub = get_info_tuple(f)
    print(io, v, " ± ", e, " ∈ [", lb, ", ", ub, "]")
end

function Base.show(io::IO, f::FitParam)
    s = Printf.@sprintf "(%.3g ± %.3g)" get_value(f) get_error(f)
    if f.frozen
        s *= "F"
    end
    print(io, s)
end

function Base.show(io::IO, ::MIME"text/plain", f::FitParam)
    print_info(io, f)
end

export FitParam,
    FitParam,
    set_value!,
    set_error!,
    get_value,
    get_error,
    get_upperlimit,
    get_lowerlimit,
    print_info
