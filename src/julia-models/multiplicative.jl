struct PhotoelectricAbsorption{T,D,F} <: AbstractTableModel{D,Multiplicative}
    table::D
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::FitParam{T}
    function PhotoelectricAbsorption(ηH::FitParam{T}) where {T}
        data = get_model_data(PhotoelectricAbsorption)
        E::Vector{Float64} = data[1]["E"]
        σ::Vector{Float64} = data[1]["σ"]
        table = LinearInterpolation(E, σ, extrapolation_bc = Line())
        new{T,typeof(table),FreeParameters{(:ηH,)}}(table, ηH)
    end
end
PhotoelectricAbsorption(; ηH = FitParam(1.0)) = PhotoelectricAbsorption(ηH)
register_model_data(PhotoelectricAbsorption, "cross_sections_phabs_angr.jld")
@fastmath function invoke!(flux, energy, ::Type{<:PhotoelectricAbsorption}, ηH, table)
    E = @views energy[1:end-1]
    @. flux = exp(-ηH * table(E))
end

export PhotoelectricAbsorption
