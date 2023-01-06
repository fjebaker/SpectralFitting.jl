
"""
    @ciskip expression

Skip the test if `CI` environment variable is test. Used to ignore tests that require downloading datasets
from the astro servers.
"""
macro ciskip(expression)
    quote
        if get(ENV, "CI", false) == false
            $expression
        else
            # skipping on CI
        end
    end |> esc
end

# add a comparison operator for models
function Base.:(==)(
    m1::SpectralFitting.AbstractSpectralModel,
    m2::SpectralFitting.AbstractSpectralModel,
)
    if typeof(m1) !== typeof(m2)
        return false
    end
    for f in fieldnames(typeof(m1))
        if getproperty(m1, f) != getproperty(m2, f)
            return false
        end
    end
    return true
end
