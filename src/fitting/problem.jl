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

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(model::FittableMultiModel))
    buff = IOBuffer()
    println(buff, "Models:")
    for m in model.m
        buf = IOBuffer()
        print(buf, "- ")
        _printinfo(buf, m)
        print(buff, indent(String(take!(buf)), 2))
    end
    print(io, encapsulate(String(take!(buff))))
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
    println(buff, "  . Models     : $(length(prob.model.m))")
    println(buff, "  . Datasets   : $(length(prob.data.d))")
    println(buff, "  Parameter Summary:")

    factor = div(length(prob.data.d), length(prob.model.m))

    total = 0
    free = 0
    frozen = 0
    for m in prob.model.m
        for p in parameter_tuple(m)
            total += factor
            if isfrozen(p)
                frozen += factor
            else
                free += factor
            end
        end
    end

    bound = 0
    for binding in prob.bindings
        bound += length(Set(first.(binding))) - 1
    end
    free = free - bound

    println(buff, "  . Total      : $(total)")
    println(
        buff,
        "  . ",
        Crayons.Crayon(foreground = :cyan),
        "Frozen",
        Crayons.Crayon(reset = true),
        "     : $(frozen)",
    )
    println(
        buff,
        "  . ",
        Crayons.Crayon(foreground = :magenta),
        "Bound",
        Crayons.Crayon(reset = true),
        "      : $(bound)",
    )
    println(
        buff,
        "  . ",
        Crayons.Crayon(foreground = :green),
        "Free (DOF)",
        Crayons.Crayon(reset = true),
        " : $(free)",
    )
    print(io, encapsulate(String(take!(buff))))
end

export FittingProblem, FittableMultiModel, FittableMultiDataset
