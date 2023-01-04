
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
