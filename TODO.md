# Converting to runtime free/frozen

Functions that need updating

    remake_with_free
    free_parameter_tuple
    free_parameter_count


## Invoke

When invoking a model, we create a `ParameterCache` to hold onto state about the FitParams. This in turn create a `Vector{T}`, which is manipulated over the course of a fit to update the free parameters


### New free

When calling a function with new free parameters
- model needs to have `FitParam` has the type
- replace element by element those which are free with new values

### New frozen

When calling a function with new frozen parameters 
- see New free

## Composite models

Number of free / frozen determined with a mix of generated and regular functions -- the generated flatten the tree and the regular functions lookup whether the parameters are free of frozen in the `FitParam`.


