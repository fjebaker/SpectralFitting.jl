push!(LOAD_PATH, "src")

using Documenter
using SpectralFitting

makedocs(
    modules = [SpectralFitting],
    clean = true,
    sitename = "SpectralFitting Documentation",
    pages = [
        "Home" => "index.md",
        "Why & How" => "why-and-how.md",
        # "Examples" => "examples.md",
        "Transitioning from XSPEC" => "transitioning-from-xspec.md",
        "Models" => "models.md",
        "Parameters" => "parameters.md",
        "Datasets" => "datasets.md",
        "Fitting" => "fitting.md"
    ],
)

deploydocs(repo = "github.com/fjebaker/SpectralFitting.jl.git")
