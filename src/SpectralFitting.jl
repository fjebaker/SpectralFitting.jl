module SpectralFitting

using LibXSPEC_jll

import Base
import Printf
import Downloads
import Statistics
import Pkg.MiniProgressBars: MiniProgressBar, start_progress, end_progress, show_progress

using FITSIO
using SparseArrays
using Surrogates
using ForwardDiff
using PreallocationTools
using MemoizedMethods
using LinearAlgebra
using FileIO
using Interpolations

import Crayons
import Parameters: @with_kw

using DocStringExtensions

include("abstract-models.jl")

include("fitparam.jl")
include("ccall-wrapper.jl")

include("composite-models.jl")
include("table-models.jl")
include("surrogate-models.jl")

include("parsing-utilities.jl")
include("function-generation.jl")

include("spectral-datasets/dataset-types.jl")
include("spectral-datasets/ogip-io.jl")
include("spectral-datasets/datasets.jl")
include("spectral-datasets/response-matrix.jl")
include("spectral-datasets/binning-utilities.jl")
include("spectral-datasets/missions/abstract-mission.jl")

include("model-data-io.jl")

include("fitting.jl")
include("plotting-recipes.jl")

# include xspec models
include("xspec-models/additive.jl")
include("xspec-models/multiplicative.jl")
include("xspec-models/convolutional.jl")

# include julia models
include("julia-models/model-utilities.jl")
include("julia-models/additive.jl")
include("julia-models/multiplicative.jl")

function __init__()
    # check if we have the minimum model data already
    _check_model_directory_present()
    # init HEASOFT
    ccall((:FNINIT, libXSFunctions), Cvoid, ())
end

end # module
