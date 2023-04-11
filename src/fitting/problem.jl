struct MultiDataset{D}
    d::D
    MultiDataset(data::Vararg{<:AbstractDataset}) = new{typeof(data)}(data)
end

domain_vector(d::MultiDataset) = reduce(vcat, (domain_vector(i) for i in d.d))
target_vector(d::MultiDataset) = reduce(vcat, (target_vector(i) for i in d.d))
target_variance(d::MultiDataset) = reduce(vcat, (target_variance(i) for i in d.d))

struct MultiModel{M}
    m::M
    MultiModel(model::Vararg{<:AbstractSpectralModel}) = new{typeof(model)}(model)
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

function _assemble_parameter_indices(bindings, n_params)
    remove = Int[]
    carry = Ref(0)
    parameter_indices = map((1:length(n_params)...,)) do i
        s = i == 1 ? 1 : n_params[i-1] + 1
        e = n_params[i]
        indices = collect(s:e)
        # assign bindings for all but the first models parameters
        if i > 1
            b = bindings[i]
            map(enumerate(indices)) do (j, index)
                # remove and replace index
                if j in b
                    push!(remove, index)
                    carry[] += 1
                    j
                else
                    # subtract the number of indices we've replaced / removed
                    index - carry[]
                end
            end
        else
            indices
        end
    end
    parameter_indices, remove
end

struct FittingProblem{M,D,B}
    model::M
    data::D
    bindings::B
    function FittingProblem(
        m::Union{<:MultiModel,AbstractSpectralModel},
        d::Union{<:MultiDataset,AbstractDataset},
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
        bindings = map(_ -> Int[], _m.m)
        new{typeof(_m),typeof(_d),typeof(bindings)}(_m, _d, bindings)
    end
end

function model_count(prob::FittingProblem)
    return length(prob.model.m)
end

function data_count(prob::FittingProblem)
    return length(prob.data.d)
end

function _get_binding_indices(prob::FittingProblem, symbols::Vararg{Symbol})
    map(prob.model.m) do model
        free_symbs =
            model isa CompositeModel ? composite_free_parameter_symbols(model) :
            free_parameter_symbols(model)
        map(symbols) do s
            i = findfirst(==(s), free_symbs)
            if isnothing(i)
                @warn "Model contains no symbol `$(Meta.quot(s))`"
                -1
            else
                i
            end
        end
    end
end

function assemble_multimodel(prob::FittingProblem)
    # unpack
    m = prob.model
    d = prob.data
    bindings = prob.bindings
    n_models = length(m.m)
    # map on tuple to keep output as tuple and therefore type stable
    funcs = map((1:n_models...,)) do i
        model = m.m[i]
        data = d.d[i]
        _lazy_folded_invokemodel(model, data)
    end
    parameters = map(freeparameters, m.m)
    all_parameters = reduce(vcat, parameters)
    impl =
        all(model -> implementation(model) isa JuliaImplementation, m.m) ?
        JuliaImplementation : XSPECImplementation

    # function which accepts all energy and parameters, and then dispatches them correctly to each sub model
    n_params = _accumulated_indices(map(length, parameters))
    n_energy = _accumulated_indices(map(data -> length(domain_vector(data)), d.d))
    n_output = _accumulated_indices(map(data -> length(target_vector(data)), d.d))

    parameter_indices, remove = _assemble_parameter_indices(bindings, n_params)
    deleteat!(all_parameters, remove)

    F = (X, params) -> begin
        fluxes = map(1:n_models) do i
            p = @views params[parameter_indices[i]]

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
        all_parameters,
        # state
        (;
            funcs = funcs,
            n_models = n_models,
            i_out = n_output,
            parameter_indices = parameter_indices,
            i_x = n_energy,
            implementation = impl,
        ),
    )
end

function bind!(prob::FittingProblem, symbols...)
    indices = _get_binding_indices(prob, symbols...)
    for (i, I) in enumerate(indices)
        if any(==(-1), I)
            j = findfirst(==(-1), I)
            throw("Could not bind symbol `$(Meta.quot(symbols[j]))`")
        end
        for index in I
            # avoid duplicating
            if index ∉ prob.bindings[i]
                push!(prob.bindings[i], index)
            end
        end
    end
    true
end

function Base.show(io::IO, ::MIME"text/plain", prob::FittingProblem)
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

export MultiDataset, MultiModel, FittingProblem, bind!
