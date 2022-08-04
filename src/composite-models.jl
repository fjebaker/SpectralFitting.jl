# model composition
struct CompositeSpectralModel{M1,M2} <: AbstractSpectralModel
    left::M1
    right::M2
end
modelkind(::Type{CompositeSpectralModel{M1,M2}}) where {M1,M2} = modelkind(M2)
invokemodel!(flux, energy, m::CompositeSpectralModel) =
    invoke_composite!(flux, energy, m::CompositeSpectralModel)

operation_string(left, right, ::Multiplicative) = "$left * $right"
operation_string(left, right, ::Convolutional) = "$left($right)"
operation_string(left, right, ::Additive) = "($left + $right)"

function modelinfo(cm::CompositeSpectralModel{M1,M2}) where {M1,M2}
    left = modelinfo(cm.left)
    right = modelinfo(cm.right)
    operation_string(left, right, modelkind(M1))
end

function Base.show(io::IO, cm::CompositeSpectralModel{M1,M2}) where {M1,M2}
    print(io, modelinfo(cm))
end
function Base.show(
    io::IO,
    ::MIME"text/plain",
    cm::CompositeSpectralModel{M1,M2},
) where {M1,M2}
    print(io, modelinfo(cm))
end

invoke_composite!(flux, energy, m::AbstractSpectralModel) = invoke!(flux, energy, m)
invoke_composite!(flux, energy, m::CompositeSpectralModel{M1,M2}) where {M1,M2} =
    invoke_composite!(flux, zeros(eltype(flux), size(flux)), energy, m)
# only used to make CompositeSpectralModel work with auto-allocations
invoke_composite!(flux, tmpflux, energy, m::AbstractSpectralModel) =
    invoke_composite!(flux, energy, m)
function invoke_composite!(flux, tmpflux, energy, m, ::Multiplicative)
    invoke_composite!(tmpflux, energy, m)
    flux .*= tmpflux
end
function invoke_composite!(flux, tmpflux, energy, m, ::Convolutional)
    invoke_composite!(flux, tmpflux, energy, m)
end
function invoke_composite!(flux, tmpflux, energy, m, ::Additive)
    invoke_composite!(tmpflux, energy, m)
    flux .+= tmpflux
end
function invoke_composite!(
    flux,
    tmpflux,
    energy,
    m::CompositeSpectralModel{M1,M2},
) where {M1,M2}
    invoke_composite!(flux, tmpflux, energy, m.right)
    invoke_composite!(flux, tmpflux, energy, m.left, modelkind(M1))
    flux
end

add_models(m1, m2, ::Additive, ::Additive) = CompositeSpectralModel(m1, m2)
add_models(_, _, ::M1, ::M2) where {M1,M2} =
    error("Left and right models must be Additive.")
add_models(m1::M1, m2::M2) where {M1,M2} = add_models(m1, m2, modelkind(M1), modelkind(M2))
function Base.:+(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel}
    add_models(m1, m2, modelkind(M1), modelkind(M2))
end

mult_models(m1, m2, ::Multiplicative, ::AbstractSpectralModelKind) =
    CompositeSpectralModel(m1, m2)
mult_models(_, _, ::M1, ::M2) where {M1,M2} = error("Left model must be Multiplicative.")
mult_models(m1::M1, m2::M2) where {M1,M2} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))
function Base.:*(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel}
    mult_models(m1, m2, modelkind(M1), modelkind(M2))
end

conv_models(m1, m2, ::Convolutional, ::Additive) = CompositeSpectralModel(m1, m2)
conv_models(_, _, ::M1, ::M2) where {M1,M2} =
    error("Left model must be Convolutional and right model must be Additive.")
conv_models(m1::M1, m2::M2) where {M1,M2} =
    conv_models(m1, m2, modelkind(M1), modelkind(M2))
function (m1::AbstractSpectralModel)(m2::M2) where {M2<:AbstractSpectralModel}
    conv_models(m1, m2)
end
