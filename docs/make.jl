push!(LOAD_PATH, "src")

using Documenter
using SpectralFitting

SpectralFitting.download_all_model_data()

makedocs(
    modules = [SpectralFitting],
    clean = true,
    sitename = "SpectralFitting.jl",
    pages = [
        "Home" => "index.md",
        "Why & How" => "why-and-how.md",
        "Examples" => "examples.md",
        # "Transitioning from XSPEC" => "transitioning-from-xspec.md",
        "Models" => [
            "Using models" => "using-models.md",
            "Model index" => "models.md",
            "Composite models" => "composite-models.md",
            "Surrogate models" => "surrogate-models.md",
        ],
        # "Parameters" => "parameters.md",
        # "Datasets" => "datasets.md",
        # "Fitting" => "fitting.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(repo = "github.com/fjebaker/SpectralFitting.jl.git")
