# Composite models

The model algebra defined by the [`AbstractSpectralModelKind`](@ref) yields instances of [`CompositeSpectralModel`](@ref), nested to various degrees. These composite models are designed to make as much information about the spectral model available at compile-time, such that rich and optimized generated functions may be assembled purely from the Julia types (see [Why & how](@ref)).

```@docs
CompositeSpectralModel
AbstractCompositeOperator
AdditionOperator
MultiplicationOperator
ConvolutionOperator
operation_symbol
```