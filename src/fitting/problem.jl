struct MultiDataset{D}
    d::D
    MultiDataset(data::Vararg{<:AbstractSpectralDataset}) = new{typeof(data)}(data)
end

struct MultiModel{M}
    m::M
    MultiModel(model::Vararg{<:AbstractSpectralModel}) = new{typeof(model)}(model)
end

struct FittingProblem{M,D}
    model::M
    data::D
    function FittingProblem(
        m::Union{<:MultiModel,AbstractSpectralModel},
        d::Union{<:MultiDataset,AbstractSpectralDataset},
    )
        _m = if !(m isa MultiModel)
            MultiModel(m)
        else
            m
        end
        _d = if !(d isa MultiDataset)
            MultiDataset(d)
        else
            d
        end
        new{typeof(_m),typeof(_d)}(_m, _d)
    end
end

function model_count(prob::FittingProblem)
    return length(prob.model.m)
end

function data_count(prob::FittingProblem)
    return length(prob.data.d)
end

function Base.show(io::IO, ::MIME"text/plain", prob::FittingProblem)
    println(io, "FittingProblem:")
    println(io, "  Models:")
    for model in prob.model.m
        println(io, " "^4 * "- $(model)")
    end
    println(io, "  Data:")
    for data in prob.data.d
        println(io, " "^4 * "- $(data)")
    end
end

export MultiDataset, MultiModel, FittingProblem
