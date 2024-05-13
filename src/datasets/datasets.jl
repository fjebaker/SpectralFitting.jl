@enumx ErrorStatistics begin
    Numeric
    Poisson
    Gaussian
    Unknown
end

abstract type AbstractDataLayout end
struct OneToOne <: AbstractDataLayout end
struct ContiguouslyBinned <: AbstractDataLayout end

export OneToOne, ContiguouslyBinned

const DEFAULT_SUPPORT_ORDERING = (ContiguouslyBinned(), OneToOne())

function preferred_support(x)
    for layout in DEFAULT_SUPPORT_ORDERING
        if supports(layout, x)
            return layout
        end
    end
    error("No prefered support for $(typeof(x))")
end

# used to check what to default to
# todo: check if can be compile time eval'd else expand for loop or @generated
function common_support(x, y)
    # order of preference is important
    # can normally trivially fallback from one-to-one to contiguous bins to regular bins
    for layout in DEFAULT_SUPPORT_ORDERING
        if supports(layout, x) && supports(layout, y)
            return layout
        end
    end
    error("No common support between $(typeof(x)) and $(typeof(y)).")
end

function _support_reducer(x::OneToOne, y)
    if supports(x, y)
        return x
    else
        error("No common support!!")
    end
end
function _support_reducer(x::ContiguouslyBinned, y)
    if supports(x, y)
        return x
    else
        _support_reducer(OneToOne(), y)
    end
end
function _support_reducer(x, y)
    common_support(x, y)
end

common_support(args::Vararg) = reduce(_support_reducer, args)

supports_contiguosly_binned(::Type) = false
supports_one_to_one(::Type) = false

supports_contiguosly_binned(x) = supports_contiguosly_binned(typeof(x))
supports_one_to_one(x) = supports_one_to_one(typeof(x))

supports(layout::AbstractDataLayout, x) = supports(layout, typeof(x))
supports(::ContiguouslyBinned, T::Type) = supports_contiguosly_binned(T)
supports(::OneToOne, T::Type) = supports_one_to_one(T)

"""
    abstract type AbstractDataset
    
Abstract type for use in fitting routines. High level representation of some underlying 
data structures. 

Fitting data is considered to have an *objective* and a *domain*. As
the domain may be, for example, energy bins (high and low), or
fourier frequencies (single value), the purpose of this abstraction
is to provide some facility for translating between these representations
for the models to fit with. This is done by checking that the [`AbstractDataLayout`](@ref)
of the model and data are compatible, or at least have compatible translations.

Must implement a minimal set of accessor methods. These are paired with
`objective` and `domain` parlance. Note that these functions are prefixed with
`make_*` and not `get_*` to represent that there may be allocations or work going
into the translation. Usage of these functions should be sparse in the interest of
performance.

The arrays returned by the `make_*` functions must correspond to the [`AbstractDataLayout`](@ref)
specified by the caller.

[`ContiguouslyBinned`](@ref) (domain should be `length(objective) + 1`, where the limits of the
``n^\text{th}`` bin are `domain[n]` and `domain[n+1]` respectively.

- [`make_objective_variance`](@ref)
- [`make_objective`](@ref)
- [`make_domain_variance`](@ref)
- [`make_model_domain`](@ref)
- [`make_ouput_domain`](@ref)

"""
abstract type AbstractDataset end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(data::AbstractDataset))
    buff = IOBuffer()
    _printinfo(buff, data)
    s = String(take!(buff))
    print(io, encapsulate(s))
end

function Base.show(io::IO, @nospecialize(data::AbstractDataset))
    print(io, "$(Base.typename(typeof(data)).name)[$(make_label(data))]")
end


"""
    make_objective

Returns the array used as the target for model fitting. The array must correspond to the data
[`AbstractDataLayout`](@ref) specified by the `layout` parameter.

In as far as it can be guarunteed, the memory in the returned array will not be mutated by any fitting procedures.

Domain for this objective should be returned by [`make_model_domain`](@ref).
"""
make_objective(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")
make_objective_variance(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

"""
    make_model_domain

Returns the array used as the domain for the modelling. This is paired with [`make_domain_variance`](@ref)
"""
make_model_domain(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")
make_domain_variance(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

"""
    make_output_domain

Returns the array used as the output domain. That is, in cases where the model
input and output map to different domains, the input domain is said to be the
model domain, the input domain is said to be the model domain. 

The distinction is mainly used for the purposes of simulating data and for
visualising data.
"""
make_output_domain(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

function objective_transformer(layout::AbstractDataLayout, dataset::AbstractDataset)
    @warn "Using default objective transformer!"
    _DEFAULT_TRANSFORMER()
end

function _DEFAULT_TRANSFORMER()
    function _transformer!!(domain, objective)
        objective
    end
    function _transformer!!(output, domain, objective)
        @. output = objective
    end
    _transformer!!
end

make_label(d::AbstractDataset) = "$(Base.typename(typeof(d)).name)"

error_statistic(::AbstractDataset) = ErrorStatistics.Unknown

"""
Must support the same API, but may also have some query methods for specific internals.
"""
abstract type AbstractMultiDataset <: AbstractDataset end

export make_model_domain,
    make_output_domain,
    make_domain_variance,
    make_objective,
    make_objective_variance,
    normalize!,
    objective_transformer

include("spectrum.jl")
include("response.jl")
include("spectraldata.jl")
include("ogipdataset.jl")
include("mission-specifics.jl")
