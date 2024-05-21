# SpectralFitting.jl Documentation

_Fast and flexible spectral fitting in Julia._

SpectralFitting.jl is a package for defining and using spectral models, with a number of utilities to make model composition easy and invocation fast. SpectralFitting wraps [LibXSPEC_jll.jl](https://github.com/astro-group-bristol/LibXSPEC_jll.jl) to expose the library of models from [HEASoft XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/), and provides helper functions for operating with spectral data from a number of different missions. The package natively uses [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) to fit parameters using the Levenberg-Marquardt algorithm, but makes it easy to use [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) for more specialized fitting algorithms, or [Turing.jl](https://github.com/TuringLang/Turing.jl) for Bayesian inference and MCMC.

SpectralFitting is designed to be extended, such that new models are simple to create, and new dataset processing pipelines for different missions are brief to define. Where performance is key, SpectralFitting helps you define fast and AD-compatible surrogates of spectral models using [Surrogates.jl](https://github.com/SciML/Surrogates.jl), and embed them in the model composition algebra.

To get started, add the AstroRegistry from the University of Bristol and then install:

```julia
julia>]
pkg> registry add https://github.com/astro-group-bristol/AstroRegistry
pkg> add SpectralFitting
```

Then use

```julia
using SpectralFitting
# ....
```

to get started. See [Walkthrough](@ref) for an example walkthrough the package.

For more University of Bristol Astrophysics Group codes, see [our GitHub organisation](https://github.com/astro-group-bristol).