module SpectralFitting

using LibXSPEC_jll

import Base
import Printf

using FITSIO
using SparseArrays
using Surrogates
using ForwardDiff
using PreallocationTools

import Crayons
import Parameters: @with_kw
import LsqFit

using DocStringExtensions

# import Turing
# import Optimization

include("abstract-models.jl")
include("fitparam.jl")
include("ccall-wrapper.jl")
include("composite-models.jl")
include("surrogate-models.jl")
# include("model-processing.jl")
# include("model-building.jl")
include("parsing-utilities.jl")
include("function-generation.jl")

include("file-io/datasets.jl")
include("file-io/response-matrix.jl")
include("file-io/binning-utilities.jl")
include("file-io/parsing-utilities.jl")
include("file-io/no-associated-mission.jl")
include("file-io/fits-loading.jl")

include("fitting.jl")
include("plotting-recipes.jl")

# include xspec models
include("xspec-models/additive.jl")
include("xspec-models/multiplicative.jl")
include("xspec-models/convolutional.jl")

# include julia models
include("julia-models/model-utilities.jl")
include("julia-models/additive.jl")

function __init__()
    # init HEASOFT
    ccall((:FNINIT, libXSFunctions), Cvoid, ())
end

end # module
