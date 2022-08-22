module SpectralUnits

using Reexport
@reexport using Unitful

Unitful.register(SpectralUnits)

@dimension 𝐄 "events" Events
@refunit counts "counts" Counts 𝐄 false

const localpromotion = Unitful.promotion
function __init__()
    Unitful.register(SpectralUnits)
    merge!(Unitful.promotion, localpromotion)
end

function infer_units(s::Symbol)
    if s == :rate
        u"counts / s"
    elseif s == :counts
        u"counts"
    else
        error("Unknown units $s.")
    end
end

export counts

const _counts = typeof(u"counts")
const _rate = typeof(u"counts / s")

end # module
