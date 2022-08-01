module SpectralFitting

using LibXSPEC_jll

using FITSIO
using DataFrames
using SparseArrays

import Crayons

import Parameters: @with_kw
import Base
import Printf

include("abstract-models.jl")
include("fit-parameters.jl")
include("ccall-wrapper.jl")
include("composite-models.jl")
include("model-processing.jl")


# include xspec models
include("xspec-models/additive.jl")
include("xspec-models/multiplicative.jl")
include("xspec-models/convolutional.jl")

function __init__()
    # init HEASOFT
    ccall((:FNINIT, libXSFunctions), Cvoid, ())
end

end # module
