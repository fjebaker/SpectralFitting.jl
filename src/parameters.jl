export FitParameter, StdFitParameter, setstd!

abstract type AbstractFitDistributionParameter <: AbstractFitParameter end

mutable struct StdFitParameter{T} <: AbstractFitParameter
    val::T
    std::T
    n_σ::Int
    frozen::Bool

    StdFitParameter(val::T, std::T; n_σ = 1, frozen::Bool = false) where {T} =
        new{T}(val, std, n_σ, frozen)
end

function Base.show(io::IO, f::StdFitParameter)
    s = Printf.@sprintf "StdFitParameter[%g ± %g" f.val f.std
    s *= (isfrozen(f) ? ", Frozen" : "") * "]"
    print(io, s)
end

upperbound(f::StdFitParameter) = f.val + f.n_σ * f.std
lowerbound(f::StdFitParameter) = f.val - f.n_σ * f.std
setupperbound!(::StdFitParameter, _) =
    error("Bound setting not supported for StdFitParameter: use `setstd!` instead.")
setlowerbound!(::StdFitParameter, _) =
    error("Bound setting not supported for StdFitParameter: use `setstd!` instead.")
setstd!(f::StdFitParameter, std) = f.std = std
std(f::StdFitParameter) = f.std

mutable struct FitParameter{T} <: AbstractFitParameter
    val::T
    lower_bound::T
    upper_bound::T
    frozen::Bool

    FitParameter(
        val::T;
        lower_bound::T = T(0.0),
        upper_bound::T = T(Inf),
        frozen::Bool = false,
    ) where {T} = new{T}(val, lower_bound, upper_bound, frozen)
end

function Base.show(io::IO, f::FitParameter)
    s = Printf.@sprintf "FitParameter[%g, lb: %g, ub: %g" f.val f.lower_bound f.upper_bound
    s *= (isfrozen(f) ? ", Frozen" : "") * "]"
    print(io, s)
end


function as_distribution(f::FitParameter)
    val = value(f)
    lb = lowerbound(f)
    ub = upperbound(f)

    (Turing.TruncatedNormal, (val, 2.0, lb, ub))
end

function as_distribution(f::StdFitParameter)
    val = value(f)
    dev = std(f)

    (Turing.Normal, (val, dev))
end
