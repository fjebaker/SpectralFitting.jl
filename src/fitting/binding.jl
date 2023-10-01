_sort_binding!(binding) = sort!(binding, by = i -> i[1])

function _construct_bound_mapping(bindings, parameter_count)
    remove = Int[]

    parameter_mapping = map((1:length(parameter_count)...,)) do i
        # indices of the full parameter array for the ith model
        collect(range(_get_range(parameter_count, i)...))
    end

    for binding in bindings
        # the model we use as a reference to bind *to*
        reference = binding[1]
        for b in @views binding[2:end]
            parameter_mapping[b[1]][b[2]] = reference[2]

            # mark for removal: find the parameter index in the global array
            N = sum(length(parameter_mapping[q]) for q = 1:b[1]-1)
            index = N + b[2]
            push!(remove, index)

            # need to now shuffle all the indices greater than this one down by 1
            for k = b[2]+1:length(parameter_mapping[b[1]])
                parameter_mapping[b[1]][k] -= 1
            end
            # and for subsequent models
            for j = b[1]+1:length(parameter_count)
                parameter_mapping[j] .-= 1
            end
        end
    end

    sort!(unique!(remove))

    parameter_mapping, remove
end

function _get_index_of_symbol(model::AbstractSpectralModel, symbol)::Int
    symbols = all_parameter_symbols(model)
    i = findfirst(==(symbol), symbols)
    if isnothing(i)
        error("Could not match symbol $symbol !")
    end
    i
end

"""

Bind the symbols of `last(pair)` in all models indexes by `first(pair)`.
"""
function _bind_pairs!(prob::FittingProblem, pairs::Vararg{<:Pair{Int,Symbol}})
    binding = map(pairs) do pair
        model = prob.model.m[pair[1]]
        pair[1] => _get_index_of_symbol(model, pair[2])
    end |> collect
    _sort_binding!(binding)
    push!(prob.bindings, binding)
end

bind!(prob::FittingProblem, pairs::Vararg{<:Pair}) = _bind_pairs!(prob, pairs...)

function bind!(prob::FittingProblem, symbs::Vararg{Symbol})
    for symb in symbs
        binding = map(enumerate(prob.model.m)) do (i, model)
            i => _get_index_of_symbol(model, symb)
        end
        push!(prob.bindings, _sort_binding!(collect(binding)))
    end
end

export bind!
