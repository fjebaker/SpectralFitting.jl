module SpectralFitting

using LibXSPEC_jll

import Base
import Printf
import Downloads
import Pkg.MiniProgressBars: MiniProgressBar, start_progress, end_progress, show_progress

using FITSIO
using SparseArrays
using Surrogates
using ForwardDiff
using PreallocationTools
using MemoizedMethods

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
include("file-io/model-data-io.jl")

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
    # check if we have the minimum model data already
    _check_model_directory_present()
    # init HEASOFT
    ccall((:FNINIT, libXSFunctions), Cvoid, ())
end

end # module
