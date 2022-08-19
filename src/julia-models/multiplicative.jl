struct PhotoelectricAbsorption{T,F1} <: AbstractTableModel{T}
    table::T
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::F1
    function PhotoelectricAbsorption(ηH::F) where {F}
        data = get_model_data(PhotoelectricAbsorption)
        E::Vector{Float64} = data[1]["E"]
        σ::Vector{Float64} = data[1]["σ"]
        table = LinearInterpolation(E, σ, extrapolation_bc=Line())
        new{typeof(table),F}(table, ηH)
    end
end
PhotoelectricAbsorption(;ηH=FitParam(1.0)) = PhotoelectricAbsorption(ηH)
register_model_data(PhotoelectricAbsorption, "cross_sections_phabs_angr.jld")
modelkind(::Type{<:PhotoelectricAbsorption}) = Multiplicative()
@fastmath function invoke!(flux, energy, ::Type{<:PhotoelectricAbsorption}, table, ηH)
    E = @views energy[1:end-1]
    @. flux = exp(-ηH * table(E))
end

export PhotoelectricAbsorption