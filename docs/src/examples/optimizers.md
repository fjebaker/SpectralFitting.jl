# Optimizer galore

Let's fit a spectrum:

```@example optimizers
using SpectralFitting

DATADIR = "..."
DATADIR = length(get(ENV, "CI", "")) > 0 ? @__DIR__() * "/../../ex-datadir" : "/home/lilith/Developer/jl/datasets/xspec/walkthrough" # hide
spec1_path = joinpath(DATADIR, "s54405.pha")
data = OGIPDataset(spec1_path) 
normalize!(data)
```