struct MultiDataset{D}
    d::D
    MultiDataset(data::Vararg{<:AbstractSpectralDataset}) = new{typeof(data)}(data)
end

energy_vector(d::MultiDataset) = reduce(vcat, (energy_vector(i) for i in d.d))
flux_vector(d::MultiDataset) = reduce(vcat, (i.rate for i in d.d))
std_dev_vector(d::MultiDataset) = reduce(vcat, (i.rateerror for i in d.d))

struct MultiModel{M}
    m::M
    MultiModel(model::Vararg{<:AbstractSpectralModel}) = new{typeof(model)}(model)
end

function _accumulated_indices(items)
    total::Int = 0
    map(items) do item
        total = item + total
        total
    end
end

function assemble_multimodel(m::MultiModel, d::MultiDataset)
    n_models = length(m.m)
    # map on tuple to keep output as tuple and therefore type stable
    funcs = map((1:n_models...,)) do i
        model = m.m[i]
        data = d.d[i]
        _lazy_folded_invokemodel(model, data)
    end
    parameters = map(freeparameters, m.m)
    autodiff =
        all(model -> implementation(model) isa JuliaImplementation, m.m) ? :forward :
        :finite

    # function which accepts all energy and parameters, and then dispatches them correctly to each sub model
    n_params = _accumulated_indices(map(length, parameters))
    n_energy = _accumulated_indices(map(data -> length(energy_vector(data)), d.d))
    n_output = _accumulated_indices(map(data -> length(data.bins_low), d.d))
    F = (X, params) -> begin
        fluxes = map(1:n_models) do i
            start_p = i == 1 ? 1 : n_params[i-1] + 1
            end_p = n_params[i]
            p = @views params[start_p:end_p]

            start_x = i == 1 ? 1 : n_energy[i-1] + 1
            end_x = n_energy[i]
            x = @views X[start_x:end_x]
            # invoke wrapped model
            funcs[i](x, p)
        end
        # combine all fluxes into single array
        reduce(vcat, fluxes)
    end

    (
        F,
        reduce(vcat, parameters),
        (;
            funcs = funcs,
            n_models = n_models,
            i_out = n_output,
            i_params = n_params,
            i_x = n_energy,
            autodiff = autodiff,
        ),
    )
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
