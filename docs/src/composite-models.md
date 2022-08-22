# Composite models

The model algebra defined by the [`AbstractSpectralModelKind`](@ref) yields instances of [`CompositeModel`](@ref), nested to various degrees. These composite models are designed to make as much information about the spectral model available at compile-time, such that rich and optimized generated functions may be assembled purely from the Julia types (see [Why & how](@ref)).

```@docs
CompositeModel
AbstractCompositeOperator
AdditionOperator
MultiplicationOperator
ConvolutionOperator
operation_symbol
```