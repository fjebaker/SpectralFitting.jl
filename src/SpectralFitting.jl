module SpectralFitting

using LibXSPEC_jll

import Base
import Printf
import Downloads
import Statistics
import Pkg.MiniProgressBars: MiniProgressBar, start_progress, end_progress, show_progress
import Distributions
import ConstructionBase

import Random

using FITSIO
using SparseArrays
using Surrogates
using LinearAlgebra
using FileIO
using Interpolations
import DataInterpolations
using SpecialFunctions
using PreallocationTools
using EnumX

import Crayons
# fitting backends
import ForwardDiff
import LsqFit
import Optimization

using DocStringExtensions

# for future use: mission specific parsing
abstract type AbstractMission end
struct NoMission <: AbstractMission end

abstract type AbstractStatistic end

# unitful units
include("units.jl")
SpectralUnits.@reexport using .SpectralUnits

include("print-utilities.jl")

include("fitparam.jl")
include("param-cache.jl")
include("abstract-models.jl")

include("ccall-wrapper.jl")

include("composite-models.jl")

include("generation/function-generation.jl")
include("generation/wrappers.jl")

include("meta-models/table-models.jl")
include("meta-models/surrogate-models.jl")

include("poisson.jl")

include("datasets/ogip.jl")
include("datasets/datasets.jl")
include("datasets/binning.jl")
include("datasets/grouping.jl")
include("datasets/injectivedata.jl")

include("model-data-io.jl")

# include fitting api
include("fitting/result.jl")
include("fitting/cache.jl")
include("fitting/problem.jl")
include("fitting/binding.jl")
include("fitting/multi-cache.jl")
include("fitting/methods.jl")
include("fitting/statistics.jl")

include("simulate.jl")
include("fitting/goodness.jl")

include("plotting-recipes.jl")

# include xspec models
include("xspec-models/additive.jl")
include("xspec-models/multiplicative.jl")
include("xspec-models/convolutional.jl")

# include julia models
include("julia-models/model-utilities.jl")
include("julia-models/additive.jl")
include("julia-models/multiplicative.jl")
include("julia-models/convolutional.jl")

function __init__()
    # check if we have the minimum model data already
    _check_model_directory_present()
    # init HEASOFT
    ccall((:FNINIT, libXSFunctions), Cvoid, ())
end

end # module
