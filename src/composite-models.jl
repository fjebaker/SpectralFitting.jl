export AbstractCompositeOperator,
    AdditionOperator,
    MultiplicationOperator,
    ConvolutionOperator,
    operation_symbol,
    CompositeSpectralModel,
    get_free_model_params,
    get_frozen_model_params

abstract type AbstractCompositeOperator end
struct AdditionOperator <: AbstractCompositeOperator end
struct MultiplicationOperator <: AbstractCompositeOperator end
struct ConvolutionOperator <: AbstractCompositeOperator end

# define these on type, since mostly used in the generation context
# which only operates on un-instantiated types
operation_symbol(::Type{<:AdditionOperator}) = :(+)
operation_symbol(::Type{<:MultiplicationOperator}) = :(*)
operation_symbol(::Type{<:ConvolutionOperator}) = nothing
operation_symbol(::O) where {O<:AbstractCompositeOperator} = operation_symbol(O)

# model composition
struct CompositeSpectralModel{M1,M2,O} <: AbstractSpectralModel
    left::M1
    right::M2
    op::O
end

modelkind(::Type{<:CompositeSpectralModel{M1,M2}}) where {M1,M2} = modelkind(M2)
function implementation(::Type{<:CompositeSpectralModel{M1,M2}}) where {M1,M2}
    if (implementation(M1) isa XSPECImplementation) ||
       (implementation(M2) isa XSPECImplementation)
        XSPECImplementation()
    else
        JuliaImplementation()
    end
end
function closurekind(::Type{<:CompositeSpectralModel{M1,M2}}) where {M1,M2}
    if (closurekind(M1) isa WithClosures) || (closurekind(M2) isa WithClosures)
        WithClosures()
    else
        WithoutClosures()
    end
end

# parameter access
# get_param_symbols(M::Type{<:AbstractSpectralModel}) = fieldnames(M)
function get_param_types(model::Type{<:SpectralFitting.CompositeSpectralModel})
    __generated_get_param_types(model)
end

# utilities
add_param_types!(_, ::Type{Nothing}) = nothing
add_param_types!(types, M::Type{<:AbstractSpectralModel}) =
    push!(types, get_param_types(M)...)

# invokation wrappers
invokemodel!(flux, energy, m::CompositeSpectralModel) =
    invoke_composite!(flux, energy, m::CompositeSpectralModel)

# composite wrappers
# only used to make CompositeSpectralModel work with auto-allocations
invoke_composite!(flux, energy, m::AbstractSpectralModel) = invokemodel!(flux, energy, m)
invoke_composite!(flux, energy, m::CompositeSpectralModel{M1,M2}) where {M1,M2} =
    invoke_composite!(flux, zeros(eltype(flux), size(flux)), energy, m)

invoke_composite!(flux, tmpflux, energy, m::AbstractSpectralModel) =
    invoke_composite!(flux, energy, m)
function invoke_composite!(flux, tmpflux, energy, m, ::Multiplicative)
    invoke_composite!(tmpflux, energy, m)
    flux .*= tmpflux
end
function invoke_composite!(flux, tmpflux, energy, m, ::Convolutional)
    invoke_composite!(flux, tmpflux, energy, m)
    flux
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


# algebra grammar
add_models(_, _, ::M1, ::M2) where {M1,M2} =
    error("Left and right models must be Additive.")
add_models(m1, m2, ::Additive, ::Additive) =
    CompositeSpectralModel(m1, m2, AdditionOperator())
add_models(m1::M1, m2::M2) where {M1,M2} = add_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:+(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel} =
    add_models(m1, m2, modelkind(M1), modelkind(M2))

mult_models(_, _, ::M1, ::M2) where {M1,M2} = error("Left model must be Multiplicative.")
mult_models(m1, m2, ::Multiplicative, ::AbstractSpectralModelKind) =
    CompositeSpectralModel(m1, m2, MultiplicationOperator())
mult_models(m1::M1, m2::M2) where {M1,M2} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:*(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))

conv_models(_, _, ::M1, ::M2) where {M1,M2} =
    error("Left model must be Convolutional and right model must be Additive.")
conv_models(m1, m2, ::Convolutional, ::Additive) =
    CompositeSpectralModel(m1, m2, ConvolutionOperator())
conv_models(m1::M1, m2::M2) where {M1,M2} =
    conv_models(m1, m2, modelkind(M1), modelkind(M2))
(m1::AbstractSpectralModel)(m2::M2) where {M2<:AbstractSpectralModel} = conv_models(m1, m2)

# runtime param access
function add_params_to_index!(params, adder!, m)
    if !isnothing(m)
        foreach(get_all_params(m)) do p
            adder!(params, p)
        end
    end
end

__get_model_parameters!(params, adder!, model::AbstractSpectralModel) =
    add_params_to_index!(params, adder!, model)

function __get_model_parameters!(params, adder!, model::CompositeSpectralModel)
    recursive_model_parse(model) do (left, right, _)
        add_params_to_index!(params, adder!, right)
        add_params_to_index!(params, adder!, left)
    end
end

function get_model_parameters(model::CompositeSpectralModel)
    params = AbstractFitParameter[]
    __get_model_parameters!(params, push!, model)
    params
end

__add_free!(ps, p, ::FreeParameter) = push!(ps, p)
__add_free!(_, _, ::FrozenParameter) = nothing
__add_free!(ps, p::F) where {F} = __add_free!(ps, p, fit_parameter_state(F))

__add_frozen!(ps, p, ::FrozenParameter) = push!(ps, p)
__add_frozen!(_, _, ::FreeParameter) = nothing
__add_frozen!(ps, p::F) where {F} = __add_frozen!(ps, p, fit_parameter_state(F))

function get_free_model_params(model::AbstractSpectralModel)
    T = generated_get_model_number_type(model)
    params = FitParam{T}[]
    __get_model_parameters!(params, __add_free!, model)
    params
end

function get_frozen_model_params(model::AbstractSpectralModel)
    T = generated_get_model_number_type(model)
    params = FrozenFitParam{T}[]
    __get_model_parameters!(params, __add_frozen!, model)
    params
end

# printing
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
