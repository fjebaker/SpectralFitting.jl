# Optimizer galore

Let's fit a spectrum:

```@example optimizers
using SpectralFitting, Plots

DATADIR = "..."
DATADIR = length(get(ENV, "CI", "")) > 0 ? @__DIR__() * "/../../ex-datadir" : "/home/lilith/Developer/jl/datasets/xspec/walkthrough" # hide
spec1_path = joinpath(DATADIR, "s54405.pha")
data = OGIPDataset(spec1_path) 
normalize!(data)

mask_energies!(data, 1, 15)

# a plotting utility
my_plot(data) = plot(
    data, 
    xscale = :log10, 
    yscale = :log10,
    ylims = (1e-3, 1.3)
)

my_plot(data)
```