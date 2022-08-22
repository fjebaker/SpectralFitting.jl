module SpectralUnits
using Reexport
@reexport using Unitful

@dimension 𝐄 "events" Events
@refunit counts "counts" Counts 𝐄 false

export counts
end
