struct FittableMultiDataset{D}
    d::D
    FittableMultiDataset(data::Vararg{<:AbstractDataset}) = new{typeof(data)}(data)
end

function Base.getindex(multidata::FittableMultiDataset, i)
    multidata.d[i]
end

struct FittableMultiModel{M}
    m::M
    FittableMultiModel(model::Vararg{<:AbstractSpectralModel}) = new{typeof(model)}(model)
end

function Base.getindex(multimodel::FittableMultiModel, i)
    multimodel.m[i]
end

function _multi_constructor_wrapper(
    T::Union{<:Type{<:FittableMultiDataset},<:Type{<:FittableMultiModel}},
    args,
)
    if args isa T
        args
    else
        if args isa Tuple
            T(args...)
        else
            T(args)
        end
    end
end

function update_model!(
    model::AbstractSpectralModel,
    result::Union{<:FittingResult,<:FittingResultSlice},
)
    for (i, f) in enumerate(_allocate_free_parameters(model))
        set_value!(f, result.u[i])
    end
    model
end

function update_model!(
    model::FittableMultiModel,
    result::Union{<:FittingResult,<:FittingResultSlice},
)
    @assert length(model.m) == 1
    update_model!(model.m[1], result)
end

function update_model!(multimodel::FittableMultiModel, result::MultiFittingResult)
    error("TODO")
end


"""
    _accumulated_indices(items)

`items` is a tuple or vector of lengths `n1, n2, ...`

Returns a tuple or array with same length as items, which gives the index boundaries of 
an array with size `n1 + n2 + ...`.
"""
function _accumulated_indices(items)
    total::Int = 0
    map(items) do item
        total = item + total
        total
    end
end

struct FittingProblem{M<:FittableMultiModel,D<:FittableMultiDataset,B}
    model::M
    data::D
    bindings::B
    function FittingProblem(m::FittableMultiModel, d::FittableMultiDataset)
        bindings = Vector{Pair{Int,Int}}[]
        new{typeof(m),typeof(d),typeof(bindings)}(m, d, bindings)
    end
end

FittingProblem(pairs::Vararg{<:Pair}) = FittingProblem(first.(pairs), last.(pairs))

function FittingProblem(m, d)
    _m = _multi_constructor_wrapper(FittableMultiModel, m)
    _d = _multi_constructor_wrapper(FittableMultiDataset, d)
    FittingProblem(_m, _d)
end

function model_count(prob::FittingProblem)
    return length(prob.model.m)
end

function data_count(prob::FittingProblem)
    return length(prob.data.d)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(prob::FittingProblem))
    buff = IOBuffer()
    println(buff, "FittingProblem:")
    println(buff, "  Models:")
    for model in prob.model.m
        println(buff, " "^4 * ". $(model)")
    end
    println(buff, "  Data:")
    for data in prob.data.d
        println(buff, " "^4 * ". $(data)")
    end
    print(io, encapsulate(String(take!(buff))))
end

export FittingProblem, FittableMultiModel, FittableMultiDataset
