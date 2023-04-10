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
using EnumX

import Crayons
# fitting backends
import ForwardDiff
import LsqFit
import Optimization

using DocStringExtensions

# unitful units
include("units.jl")
include("print-utilities.jl")

include("fitparam.jl")
include("abstract-models.jl")

include("ccall-wrapper.jl")

include("composite-models.jl")

include("generation/function-generation.jl")
include("generation/wrappers.jl")

include("meta-models/table-models.jl")
# include("meta-models/surrogate-models.jl")

include("poisson.jl")

include("datasets/types.jl")
include("datasets/ogip.jl")
include("datasets/spectral-datasets.jl")
include("datasets/response-matrix.jl")
include("datasets/binning-utilities.jl")
include("datasets/missions/abstract-mission.jl")
include("datasets/simple-dataset.jl")

include("model-data-io.jl")


# include fitting api
include("fitting/problem.jl")
include("fitting/result.jl")
include("fitting/methods.jl")
include("fitting/statistics.jl")

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
