"""
    SurrogateSpectralModel <: AbstractSpectralModel
    SurrogateSpectralModel(modelkind, surrogate, params, params_symbols)

Used to wrap a surrogate function into an [`AbstractSpectralModel`](@ref).

# Example

Creating a surrogate function using [`make_surrogate_function`](@ref):
```julia
# build and optimize a surrogate model
surrogate = make_surrogate_function(model, lower_bounds, upper_bounds)

# create surrogate spectral model
sm = SurrogateSpectralModel(
    Multiplicative(),
    surrogate,
    (FitParam(1.0),),
    (:ηH,)
)
```

The `lower_bounds` and `upper_bounds` must be tuples in the form `(E, params...)`, where `E`
denotes the bounds on the energy range to train over.
"""
struct SurrogateSpectralModel{T,K,N,S,Symbols} <: AbstractSpectralModel{T,K}
    surrogate::S
    params::NTuple{N,T}
end

function SurrogateSpectralModel(
    ::K,
    surrogate::S,
    params::Tuple,
    symbols,
) where {K<:AbstractSpectralModelKind,S}
    P = typeof(params)
    N = length(P.parameters)
    T = first(P.parameters)
    SurrogateSpectralModel{T,K,N,S,symbols}(surrogate, params)
end


# reflection tie-ins
function Reflection.get_closure_symbols(::Type{<:SurrogateSpectralModel})
    (:surrogate,)
end
function Reflection.get_parameter_symbols(M::Type{<:SurrogateSpectralModel})
    last(M.parameters)
end
function Reflection.parameter_lenses(
    ::Type{<:SurrogateSpectralModel},
    info::Reflection.ModelInfo,
)
    map(eachindex(info.symbols)) do i
        :(getfield($(info.lens), :params)[$i])
    end
end
function Reflection.make_constructor(
    M::Type{<:SurrogateSpectralModel},
    closures::Vector,
    params::Vector,
    T::Type,
)
    _, K, N, S, Syms = M.parameters
    :(SurrogateSpectralModel{$T,$K,$N,$S,$Syms}($(closures...), ($(params...),)))
end

function remake_with_number_type(
    model::SurrogateSpectralModel{FitParam{T},K,N,S,Syms},
) where {T,K,N,S,Syms}
    params = model_parameters_tuple(model)
    new_params = convert.(T, params)
    SurrogateSpectralModel{T,K,N,S,Syms}(model.surrogate, NTuple{N,T}(new_params))
end

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

@fastmath function invoke!(
    flux,
    energy,
    model::SurrogateSpectralModel{T,<:Multiplicative},
) where {T}
    @inbounds for i in eachindex(flux)
        E = T(energy[i])
        v = (E, model.params...)
        flux[i] = model.surrogate(v)
    end
end

@fastmath function invoke!(
    flux,
    energy,
    model::SurrogateSpectralModel{T,<:Additive},
) where {T}
    finite_diff_kernel!(flux, energy) do E
        model.surrogate((E, model.params...))
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
function wrap_model_as_objective(model::AbstractSpectralModel; ΔE = 1e-1)
    (x) -> begin
        energies = [first(x), first(x) + ΔE]
        flux = zeros(typeof(x[2]), 1)
        invokemodel!(flux, energies, model, [x[2:end]...]) |> first
    end
end

function _initial_space(obj, lb, ub, sample_type, N)
    xys = sample(N, lb, ub, sample_type)
    zs = obj.(xys)
    (; x = xys, y = zs)
end

struct SurrogateHarness{M,O,S,L}
    model::M
    obj::O
    surrogate::S
    lower_bounds::L
    upper_bounds::L
end

function optimize_accuracy!(harness::SurrogateHarness, ; kwargs...)
    optimize_accuracy!(
        harness.surrogate,
        harness.obj,
        harness.lower_bounds,
        harness.upper_bounds;
        kwargs...,
    )
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

!!! warning
    Additive models integrate energies to calculate flux, which surrogate models are currently not
    capable of. Results for Additive models likely to be inaccurate. This will be patched in a future
    version.
"""
function make_surrogate_function(
    make_surrogate::Function,
    model::M,
    lowerbounds::T,
    upperbounds::T,
    ;
    seed_samples = 50,
    sample_type = SobolSample(),
) where {M<:AbstractSpectralModel,T<:NTuple}
    # do tests here to make sure lower and upper bounds are okay
    obj = wrap_model_as_objective(model)

    K = modelkind(M)
    if K === Additive()
        @warn "Additive models integrate energies to calculate flux, which spectral surrogate models currently do not do. Surrogate is consequently likely to have poor accuracy for Additive models."
    end

    # build surrogate
    inits = _initial_space(obj, lowerbounds, upperbounds, sample_type, seed_samples)

    @show typeof(inits)

    SurrogateHarness(model, obj, make_surrogate(inits.x, inits.y), lowerbounds, upperbounds)
end

function make_model(harness::SurrogateHarness)
    SurrogateSpectralModel(
        modelkind(harness.model),
        harness.surrogate,
        model_parameters_tuple(harness.model),
        all_parameter_symbols(harness.model),
    )
end

export SurrogateSpectralModel,
    optimize_accuracy!, wrap_model_as_objective, make_surrogate_function, make_model
