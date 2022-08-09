"""
    SurrogateSpectralModel <: AbstractSpectralModel
    SurrogateSpectralModel(modelkind, surrogate, params, params_symbols)

Used to wrap a surrogate function into an [`AbstractSpectralModel`](@ref).

# Example

Creating a surrogate function for [`XS_PhotoelectricAbsorption`](@ref) using
[`wrap_model_as_objective`](@ref) and [`make_surrogate_function`](@ref).
```julia
using SpectralFitting
using Surrogates

model = XS_PhotoelectricAbsorption()
objsampler = wrap_model_as_objective(model)

# upper and lower bounds on energy and ηH
lb = (0.1, 1e-3)
ub = (30.0, 30.0)

# build surrogate function
surrogate = make_surrogate_function(model, s_lb, s_ub)

# create surrogate model
sm = SurrogateSpectralModel(
    Multiplicative(),
    surrogate,
    (FitParam(1.0),),
    (:ηH,)
)
```
"""
struct SurrogateSpectralModel{K,S,P,Z} <: AbstractSpectralModel
    surrogate::S
    params::P
    params_symbols::Z
    SurrogateSpectralModel(
        ::K,
        surrogate::S,
        params::P,
        params_symbols::Z,
    ) where {K,S,P,Z} = new{K,S,P,Z}(surrogate, params, params_symbols)
end

closurekind(::Type{<:SurrogateSpectralModel}) = WithClosures()
model_base_name(::Type{<:SurrogateSpectralModel{K}}) where {K} =
    :(SurrogateSpectralModel{$K})

# model generation
get_closure_param_fields(::Type{<:SurrogateSpectralModel}) = (:surrogate,)
get_param_types(::Type{<:SurrogateSpectralModel{K,S,P,Z}}) where {K,S,P,Z} = P.types
get_param_symbols(M::Type{<:SurrogateSpectralModel}) =
    [Symbol(:P, i) for i in eachindex(get_param_types(M))]

# runtime access
get_param_symbols(m::SurrogateSpectralModel) = m.params_symbols
function get_param(m::SurrogateSpectralModel, s::Symbol)
    i = findfirst(==(s), m.params_symbols)
    if !isnothing(i)
        m.params[i]
    else
        error("No such symbol: $s")
    end
end

modelkind(::Type{<:SurrogateSpectralModel{Additive}}) = Additive()
modelkind(::Type{<:SurrogateSpectralModel{Multiplicative}}) = Multiplicative()
modelkind(::Type{<:SurrogateSpectralModel{Convolutional}}) = Convolutional()

@fastmath function invoke!(
    flux,
    energy,
    ::Type{<:SurrogateSpectralModel{Multiplicative}},
    surrogate,
    params::T...,
) where {T}
    @inbounds for i in eachindex(flux)
        E = T(energy[i])
        v = (E, params...)
        flux[i] = surrogate(v)
    end
end

@fastmath function invoke!(
    flux,
    energy,
    ::Type{<:SurrogateSpectralModel{Additive}},
    surrogate,
    params...,
)
    integrate_over_flux!(flux, energy) do E
        surrogate((E, params...))
    end
end

"""
    optimize_accuracy!(
        surr::AbstractSurrogate,
        obj::Function,
        lb,
        ub;
        sample_type::SamplingAlgorithm = SobolSample(),
        maxiters = 200,
        N_truth = 5000,
        verbose = false,
    )

Improve accuracy (faithfullness) of the surrogate model in recreating the objective function.

Samples a new space of `N_truth` points between `lb` and `ub`, and calculates the objective
function `obj` at each. Finds the point with largest MSE between surrogate and objective, and
adds the point to the surrogate pool. Repeats `maxiters` times, adding `maxiters` points to
surrogate model.

Optionally print to stdout the MSE and iteration count with `verbose = true`.

Note that upper- and lower-bounds should be in the form `(E, params...)`, where `E` is a single
energy and `params` are the model parameters.
"""
function optimize_accuracy!(
    surr::AbstractSurrogate,
    obj::Function,
    lb,
    ub;
    sample_type::SamplingAlgorithm = SobolSample(),
    maxiters = 200,
    N_truth = 5000,
    verbose = false,
)
    X = sample(N_truth, lb, ub, sample_type)
    y = obj.(X)
    for epoch = 1:maxiters
        ŷ = surr.(X)
        σ = @. (ŷ - y)^2
        ℳσ, i = findmax(σ)
        if verbose
            println("$(rpad(epoch, 5)): ", ℳσ)
        end
        add_point!(surr, X[i], y[i])
    end
    surr
end

"""
    wrap_model_as_objective(model::AbstractSpectralModel; ΔE = 1e-1)
    wrap_model_as_objective(M::Type{<:AbstractSpectralModel}; ΔE = 1e-1)

Wrap a spectral model into an objective function for building/optimizing a surrogate model.
Returns an anonymous function taking the tuple `(E, params...)` as the argument, and
returning a single flux value.
"""
wrap_model_as_objective(::M; kwargs...) where {M<:AbstractSpectralModel} =
    wrap_model_as_objective(M; kwargs...)
function wrap_model_as_objective(M::Type{<:AbstractSpectralModel}; ΔE = 1e-1)
    (x) -> begin
        energies = [first(x), first(x) + ΔE]
        flux = zero(typeof(x[2]))
        invokemodel!(flux, energies, M, x[2:end]...)[1]
    end
end
function wrap_model_as_objective(model::CompositeSpectralModel; ΔE = 1e-1)
    (x) -> begin
        energies = [first(x), first(x) + ΔE]
        flux = make_fluxes(energies, flux_count(model), typeof(x[2]))
        generated_model_call!(flux, energies, model, x[2:end])[1]
    end
end

function _initial_space(obj, lb, ub, sample_type, N)
    xys = sample(N, lb, ub, sample_type)
    zs = obj.(xys)
    (xys, zs)
end

"""
    make_surrogate_function(
        model::M,
        lowerbounds::T,
        upperbounds::T;
        optimization_samples = 200,
        seed_samples = 50,
        S::Type = RadialBasis,
        sample_type = SobolSample(),
        verbose = false,
    )

Creates and optimizes a surrogate model of type `S` for `model`, using [`wrap_model_as_objective`](@ref)
and  [`optimize_accuracy!`](@ref) for `optimization_samples` iterations. Model is initially
seeded with `seed_samples` points prior to optimization.

!!! note
    Additive models integrate energies to calculate flux, which surrogate models are currently not
    capable of. Results for Additive models likely to be inaccurate. This will be patched in a future
    version.
"""
function make_surrogate_function(
    model::M,
    lowerbounds::T,
    upperbounds::T;
    optimization_samples = 200,
    seed_samples = 50,
    S::Type = RadialBasis,
    sample_type = SobolSample(),
    verbose = false,
) where {M<:AbstractSpectralModel,T<:NTuple{2,N}} where {N}
    # do tests here to make sure lower and upper bounds are okay
    obj = wrap_model_as_objective(model)

    K = modelkind(M)
    if K === Additive()
        @warn "Additive models integrate energies to calculate flux, which spectral surrogate models currently do not do. Surrogate is consequently likely to have poor accuracy for Additive models."
    end

    # build surrogate
    (X::Vector{T}, Y::Vector{N}) =
        _initial_space(obj, lowerbounds, upperbounds, sample_type, seed_samples)
    surrogate = S(X, Y, lowerbounds, upperbounds)
    # optimize surrogate
    optimize_accuracy!(
        surrogate,
        obj,
        lowerbounds,
        upperbounds;
        sample_type = sample_type,
        maxiters = optimization_samples,
        verbose = verbose,
    )
    surrogate
end

export SurrogateSpectralModel,
    optimize_accuracy!, wrap_model_as_objective, make_surrogate_function
