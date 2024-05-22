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
        frozen = false,
        lower_limit = T(0.0),
        upper_limit = T(Inf),
        error = fit_param_default_error(val),
    ) where {T} = new{T}(val, error, lower_limit, upper_limit, frozen)
end

# interface
set_value!(f::FitParam{T}, val::T) where {T} = f.value = val
set_error!(f::FitParam{T}, val::T) where {T} = f.error = val
get_value(f::FitParam) = f.value
isfrozen(f::FitParam) = f.frozen
isfree(f::FitParam) = !isfrozen(f)
function set!(f::FitParam, o::FitParam)
    f.value = o.value
    f.lower_limit = o.lower_limit
    f.upper_limit = o.upper_limit
    f.error = o.error
    f
end
# edge case
get_value(x::Number) = x

get_error(f::FitParam) = f.error
get_upperlimit(f::FitParam) = f.upper_limit
get_lowerlimit(f::FitParam) = f.lower_limit

Base.isapprox(f1::FitParam, f2::FitParam; kwargs...) =
    isapprox(f1.value, f2.value; kwargs...)
Base.:(==)(f1::FitParam, f2::FitParam) = f1.value == f2.value
Base.convert(T::Type{<:Number}, f::FitParam) = convert(T, f.value)

paramtype(::Type{FitParam{T}}) where {T} = T
paramtype(::T) where {T<:FitParam} = paramtype(T)

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

function Base.show(io::IO, @nospecialize(f::FitParam))
    s = Printf.@sprintf "(%.3g ± %.3g)" get_value(f) get_error(f)
    print(io, s)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(f::FitParam))
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
