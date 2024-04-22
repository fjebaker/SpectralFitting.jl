_sort_binding!(binding) = sort!(binding, by = i -> i[1])

function _construct_bound_mapping(bindings, parameter_count)
    println("Bindings are ", bindings)
    println("Parameter count is ", parameter_count, " with length ", length(parameter_count))
    remove = Int[]

    parameter_mapping = map((1:length(parameter_count)...,)) do i
        # indices of the full parameter array for the ith model
        collect(range(_get_range(parameter_count, i)...))
    end
    println("Initial parameter mapping is ", parameter_mapping)

    for binding in bindings
        println("  Looking at binding ", binding)
        # the model we use as a reference to bind *to*
        reference = binding[1]
        println("  Reference to bind to is ", reference)
        for b in @views binding[2:end]
            println("    Now considering binding ", b, " which has elements ", b[1], " and ", b[2])
            # get the number of the parameter that we are binding
            parameter_number = parameter_mapping[b[1]][b[2]]
            parameter_mapping[b[1]][b[2]] = reference[2]
            println("    Reference[2] is ", reference[2])
            println("    Parameter mapping is ", parameter_mapping)

            # mark for removal: find the parameter index in the global array
            N = sum(length(parameter_mapping[q]) for q = 1:b[1]-1)
            index = N + b[2]
            push!(remove, index)
            println("    Remove updated to be ", remove, " having added ", index, " with value ", parameter_number)

            # need to now shuffle all the indices greater than this one down by 1
            for k = b[2]+1:length(parameter_mapping[b[1]])
                if (parameter_mapping[b[1]][k] > parameter_number)
                    parameter_mapping[b[1]][k] -= 1
                    println("      Changed parameter_mapping[", b[1], "][", k, "] to ", parameter_mapping[b[1]][k])
                end
            end
            # and for subsequent models
            for j = b[1]+1:length(parameter_count)
                for k = 1:length(parameter_mapping[j])
                    if parameter_mapping[j][k] > parameter_number
                        parameter_mapping[j][k] -= 1
                        println("      Also changed mapping[", j, "][", k, "] to ", parameter_mapping[j][k])
                    end
                end
            end
            println("    Parameter mapping is ", parameter_mapping)
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
