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
            # get the number of the parameter that we are binding
            parameter_number = parameter_mapping[b[1]][b[2]]
            parameter_mapping[b[1]][b[2]] = reference[2]

            # mark for removal: find the parameter index in the global array
            N = if b[1] > 1
                sum(length(parameter_mapping[q]) for q = 1:b[1]-1)
            else
                0
            end
            index = N + b[2]
            push!(remove, index)

            # need to now shuffle all the indices greater than this one down by 1
            for k = b[2]+1:length(parameter_mapping[b[1]])
                if (parameter_mapping[b[1]][k] > parameter_number)
                    parameter_mapping[b[1]][k] -= 1
                end
            end
            # and for subsequent models
            for j = b[1]+1:length(parameter_count)
                for k = 1:length(parameter_mapping[j])
                    if parameter_mapping[j][k] > parameter_number
                        parameter_mapping[j][k] -= 1
                    end
                end
            end
        end
    end

    sort!(unique!(remove))

    parameter_mapping, remove
end

function _get_index_of_symbol(model::AbstractSpectralModel, symbol)::Int
    pnt = parameter_named_tuple(model)
    symbols = keys(pnt)
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

"""
    bind!(prob::FittingProblem, pairs::Pair{Int,Symbol}...)
    bind!(prob::FittingProblem, symbols::Symbol...)

Bind parameters together within a [`FittingProblem`](@ref). Parameters bound
together will be mandated to have exact same value during the fit.

The binding may either be a single symbol that is present in all models in the
fitting problem, or a series of pairs `Int => Symbol` which index the specific
model and parameters to bind together. All bindings specified in a single call
to `bind!` will be bound together. Multiple bindings are possible with repeated
call to `bind!`.

- Bind model 1's `K_1` parameter to model 2's `K_3`:

  ```julia
  bind!(prob, 1 => :K_1, 2 => :K_3)
  ```

- Bind model 3's `K_2` parameter to model4's `:L_1` and model 6's `a_3`:

  ```julia
  bind!(prob, 3 => :K_2, 4 => :L_1, 6 => :a_3)
  ```

- Bind the `K_1` parameter across all the models:

  ```julia
  bind!(prob, :K_1)
  ```

## Examples

Consider the following two models
```julia
model1 = PhotoelectricAbsorption() * (BlackBody() + PowerLaw())
model2 = PhotoelectricAbsorption() * (PowerLaw() + PowerLaw())

prob = FittingProblem(model1 => data1, model2 => data2)

# bind the power law indices in the two models
bind!(prob, :a_1)

# bind the normalisation of powerlaws in the 2nd model:
bind!(prob, 2 => :K_1, 2 => :K_2)
```

!!! note
    Only free parameters can be bound together.
"""
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
