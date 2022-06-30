module SpectralFitting

using LibXSPEC_jll
using Parameters

include("abstract-models.jl")
include("ccall-wrapping.jl")

for xspec_model in readdir(joinpath(@__DIR__, "xspec-models"); join = true)
    if last(splitext(xspec_model)) == ".jl"
        include(xspec_model)
    end
end

for model in readdir(joinpath(@__DIR__, "julia-models"); join = true)
    if last(splitext(model)) == ".jl"
        include(model)
    end
end

function __init__()
    # init HEASOFT
    ccall((:FNINIT, libXSFunctions), Cvoid, ())
end

end # module
