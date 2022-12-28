"""
Utility structure for tracking which parameters are free or frozen.
"""
struct FreeParameters{V} end

function fit_param_default_error(val)
    # 10 % error
    round(abs(0.1 * val), sigdigits = 1)
end

abstract type AbstractParameterState end
struct FrozenParameter <: AbstractParameterState end
struct FreeParameter <: AbstractParameterState end

is_frozen(::FrozenParameter) = true
is_frozen(::FreeParameter) = false

abstract type AbstractFitParameter end
is_frozen(::F) where {F<:AbstractFitParameter} = is_frozen(F)
is_frozen(F::Type{<:AbstractFitParameter}) = is_frozen(fit_parameter_state(F))

# interface

set_value!(f::AbstractFitParameter, val) = f.value = val
set_error!(f::AbstractFitParameter, val) = f.error = val

get_value(x::Number) = x
get_value(f::AbstractFitParameter) = f.value

get_error(f::AbstractFitParameter) = f.error
get_upperlimit(f::AbstractFitParameter) = f.upper_limit
get_lowerlimit(f::AbstractFitParameter) = f.lower_limit

fit_parameter_state(::Type{<:AbstractFitParameter}) = FreeParameter()
fit_parameter_state(::F) where {F<:AbstractFitParameter} = fit_parameter_state(F)

as_distribution(f::AbstractFitParameter) =
    (Turing.TruncatedNormal, (f.value, f.error, f.lower_limit, f.upper_limit))

# concrete types

mutable struct FitParam{T} <: AbstractFitParameter
    value::T
    error::T

    lower_limit::T
    upper_limit::T

    FitParam(
        val::T;
        lower_limit = T(0.0),
        upper_limit = T(Inf),
        error = fit_param_default_error(val),
    ) where {T} = new{T}(val, error, lower_limit, upper_limit)
end

Base.isapprox(f1::FitParam, f2::FitParam; kwargs...) = isapprox(f1.value, f2.value; kwargs...)

parameter_type(::Type{FitParam{T}}) where {T} = T
parameter_type(::T) where {T<:FitParam} = parameter_type(T)

function get_info_tuple(f::AbstractFitParameter)
    s1 = Printf.@sprintf "%.3g" get_value(f)
    s2 = Printf.@sprintf "%.3g" get_error(f)
    s3 = Printf.@sprintf "%.3g" get_lowerlimit(f)
    s4 = Printf.@sprintf "%.3g" get_upperlimit(f)
    (s1, s2, s3, s4)
end

function print_info(io::IO, f::AbstractFitParameter)
    v, e, lb, ub = get_info_tuple(f)
    f = if is_frozen(f)
        Crayons.Box.CYAN_FG("Frozen")
    else
        Crayons.Box.GREEN_FG("Free")
    end
    print(io, v, " ± ", e, " ∈ [", lb, ", ", ub, "] ", f)
end

function Base.show(io::IO, f::AbstractFitParameter)
    s = Printf.@sprintf "(%.3g ± %.3g)" get_value(f) get_error(f)
    print(io, s)
end

function Base.show(io::IO, ::MIME"text/plain", f::AbstractFitParameter)
    print_info(io, f)
end

export AbstractFitParameter,
    FitParam,
    set_value!,
    set_error!,
    get_value,
    get_error,
    get_upperlimit,
    get_lowerlimit,
    is_frozen,
    as_distribution,
    print_info
