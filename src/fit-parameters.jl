mutable struct FitParameter{T} <: AbstractFitParameter
    val::T
    lower_bound::T
    upper_bound::T
    frozen::Bool

    FitParameter(val::T; lower_bound::T = T(0.0), upper_bound::T = T(Inf), frozen::Bool = false) where {T} = new{T}(val, lower_bound, upper_bound, frozen)
end

function Base.show(io::IO, f::FitParameter)
    print(io,
        "FitParam[$(f.val), lb: $(f.lower_bound), ub: $(f.upper_bound)" * (f.frozen ? ", Frozen" : "")  * "]"
    )
end
