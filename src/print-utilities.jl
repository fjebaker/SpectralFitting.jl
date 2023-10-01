
function prettyfloat(f)
    if f == 0
        "0.0"
    elseif f == Inf
        "Inf"
    elseif (f ≥ 1) && (f - trunc(Int, f) < 1e-8)
        Printf.@sprintf("%.1f", f)
    else
        Printf.@sprintf("%#.5g", f)
    end
end

function encapsulate(text)
    # drop everything after last new line
    n = findlast(==('\n'), text)
    s = if strip(text[n:end]) == ""
        text[1:n]
    else
        text
    end
    lines = split(s, "\n")
    out = map(enumerate(lines)) do (i, line)
        if i == 1
            "┌ " * line
        elseif i == lastindex(lines)
            "└ " * line
        else
            "│ " * line
        end
    end
    join(out, "\n")
end

function indent(text, n)
    spacer = " "^n
    replace(text, '\n' => "\n" * spacer)
end
