# model composition
struct CompositeSpectralModel{M1,M2} <: AbstractSpectralModel
    left::M1
    right::M2
end
modelkind(::Type{CompositeSpectralModel{M1,M2}}) where {M1,M2} = modelkind(M2)
invokemodel!(flux, energy, m::CompositeSpectralModel) =
    __invoke_composite!(flux, energy, m::CompositeSpectralModel)

function modelinfo(cm::CompositeSpectralModel{M1,M2}) where {M1,M2}
    left = modelinfo(cm.left)
    right = modelinfo(cm.right)
    if modelkind(M1) === Multiplicative
        "$left * $right"
    elseif modelkind(M1) === Additive
        "($left + $right)"
    else
        "$left($right)"
    end
end

function Base.show(io::IO, cm::CompositeSpectralModel{M1,M2}) where {M1,M2}
    print(io, modelinfo(cm))
end

__invoke_composite!(flux, energy, m::AbstractSpectralModel) = invoke!(flux, energy, m)
__invoke_composite!(flux, energy, m::CompositeSpectralModel{M1,M2}) where {M1,M2} =
    __invoke_composite!(flux, zeros(eltype(flux), size(flux)), energy, m)

# only used to make CompositeSpectralModel work with auto-allocations
__invoke_composite!(flux, tmpflux, energy, m::AbstractSpectralModel) =
    __invoke_composite!(flux, energy, m)
# this could be a generated function
function __invoke_composite!(
    flux,
    tmpflux,
    energy,
    m::CompositeSpectralModel{M1,M2},
) where {M1,M2}
    __invoke_composite!(flux, tmpflux, energy, m.right)

    if modelkind(M1) === Convolutional
        __invoke_composite!(flux, tmpflux, energy, m.left)

    elseif modelkind(M1) === Multiplicative
        __invoke_composite!(tmpflux, energy, m.left)
        flux .*= tmpflux

    elseif modelkind(M1) === Additive
        __invoke_composite!(tmpflux, energy, m.left)
        flux .+= tmpflux

    else
        error("Unknown model kind.")

    end

    flux
end

function Base.:+(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel}
    if (modelkind(M1) === Additive) && (modelkind(M2) === Additive)
        CompositeSpectralModel(m1, m2)
    else
        error("Left and right models must be Additive.")
    end
end

function Base.:*(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel}
    if (modelkind(M1) === Multiplicative) && (modelkind(M2) === Additive)
        CompositeSpectralModel(m1, m2)
    else
        error("Left model must be Multiplicative and right model must be Additive.")
    end
end

function (m1::AbstractSpectralModel)(m2::M2) where {M2<:AbstractSpectralModel}
    if (modelkind(typeof(m1)) === Convolutional) && (modelkind(M2) === Additive)
        CompositeSpectralModel(m1, m2)
    else
        error("Left model must be Convolutional and right model must be Additive.")
    end
end
