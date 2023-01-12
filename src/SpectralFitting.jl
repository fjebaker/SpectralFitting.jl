module SpectralFitting

using LibXSPEC_jll

import Base
import Printf
import Downloads
import Statistics
import Pkg.MiniProgressBars: MiniProgressBar, start_progress, end_progress, show_progress
import Distributions
import ConstructionBase

using FITSIO
using SparseArrays
using Surrogates
using LinearAlgebra
using FileIO
using Interpolations
using SpecialFunctions

import Crayons
# fitting backends
import LsqFit
import Optimization
import OptimizationOptimJL

using DocStringExtensions

# unitful units
include("units.jl")

include("fitparam.jl")
include("abstract-models.jl")

include("ccall-wrapper.jl")

include("composite-models.jl")

include("generation/function-generation.jl")
include("generation/wrappers.jl")

include("meta-models/table-models.jl")
# include("meta-models/surrogate-models.jl")

include("poisson.jl")

include("spectral-datasets/dataset-types.jl")
include("spectral-datasets/ogip-io.jl")
include("spectral-datasets/datasets.jl")
include("spectral-datasets/response-matrix.jl")
include("spectral-datasets/binning-utilities.jl")
include("spectral-datasets/missions/abstract-mission.jl")

include("model-data-io.jl")


# include fitting api
include("fitting/problem.jl")
include("fitting/result.jl")
include("fitting/methods.jl")

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
