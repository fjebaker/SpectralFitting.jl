struct PhotoelectricAbsorption{D,T,F} <: AbstractTableModel{T,Multiplicative}
    table::D
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::T
end
function PhotoelectricAbsorption(ηH::T) where {T}
    data = get_model_data(PhotoelectricAbsorption)
    E::Vector{Float64} = data[1]["E"]
    σ::Vector{Float64} = data[1]["σ"]
    table = LinearInterpolation(E, σ, extrapolation_bc = Line())
    PhotoelectricAbsorption{typeof(table),T,FreeParameters{(:ηH,)}}(table, ηH)
end
PhotoelectricAbsorption(; ηH = FitParam(1.0)) = PhotoelectricAbsorption(ηH)
register_model_data(PhotoelectricAbsorption, "cross_sections_phabs_angr.jld")
@inline @fastmath function invoke!(flux, energy, model::PhotoelectricAbsorption)
    let ηH = model.ηH, table = model.table
        E = @views energy[1:end-1]
        @. flux = exp(-ηH * table(E))
    end
end

export PhotoelectricAbsorption
