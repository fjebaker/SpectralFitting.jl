function fit_param_default_error(val)
    # 10 % error
    round(abs(0.1 * val), sigdigits = 1)
end

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

function get_info_tuple(f::FitParam)
    err = get_error(f)
    s1 = if abs(err) > eps(typeof(err))
        Printf.@sprintf "%.3g ± %.3g" get_value(f) err
    else
        Printf.@sprintf "%.3g" get_value(f)
    end
    s2 = Printf.@sprintf "(%.3g - %.3g)" get_lowerlimit(f) get_upperlimit(f)
    (s1, s2)
end

function print_info(io::IO, f::FitParam)
    s1, s2 = get_info_tuple(f)
    s3 = if is_frozen(f)
        Crayons.Box.CYAN_FG("Frozen")
    else
        Crayons.Box.GREEN_FG("Free")
    end
    print(io, s1, " ∈ ", s2, " ", s3)
end

function Base.show(io::IO, f::FitParam)
    s = Printf.@sprintf "(%.3g ± %.3g)" get_value(f) get_error(f)
    print(io, s)
end

function Base.show(io::IO, ::MIME"text/plain", f::FitParam)
    print_info(io, f)
end

# interface

set_value!(f::FitParam, val) = f.value = val
set_error!(f::FitParam, val) = f.error = val

get_value(x::Number) = x
get_value(f::FitParam) = f.value

get_error(f::FitParam) = f.error
get_upperlimit(f::FitParam) = f.upper_limit
get_lowerlimit(f::FitParam) = f.lower_limit


is_frozen(f::FitParam) = f.frozen
set_freeze!(f::FitParam, state) = f.frozen = state
freeze!(f::FitParam) = set_freeze!(f, true)
unfreeze!(f::FitParam) = set_freeze!(f, false)

as_distribution(f::FitParam) =
    (Turing.TruncatedNormal, (f.value, f.error, f.lower_limit, f.upper_limit))


export FitParam,
    set_value!,
    set_error!,
    get_value,
    get_error,
    get_upperlimit,
    get_lowerlimit,
    is_frozen,
    set_freeze!,
    freeze!,
    unfreeze!,
    as_distribution,
    print_info
