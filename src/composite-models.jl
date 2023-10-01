export AbstractCompositeOperator,
    AdditionOperator,
    MultiplicationOperator,
    ConvolutionOperator,
    operation_symbol,
    CompositeModel

"""
    abstract type AbstractCompositeOperator

Superype of all composition operators. Used to implement the model algebra of [`AbstractSpectralModelKind`](@ref)
through a trait system.

The Julia symbol corresponding to a given `AbstractCompositeOperator` may be obtained through
[`operation_symbol`](@ref).
"""
abstract type AbstractCompositeOperator end
"""
    AdditionOperator <: AbstractCompositeOperator
    AdditionOperator()

Corresponds to the `:(+)` symbol.
"""
struct AdditionOperator <: AbstractCompositeOperator end
"""
    MultiplicationOperator <: AbstractCompositeOperator
    MultiplicationOperator()

Corresponds to the `:(*)` symbol.
"""
struct MultiplicationOperator <: AbstractCompositeOperator end
"""
    ConvolutionOperator <: AbstractCompositeOperator
    ConvolutionOperator()

Has no corresponding symbol, since it invokes a function call `C(A)`.
"""
struct ConvolutionOperator <: AbstractCompositeOperator end

# define these on type, since mostly used in the generation context
# which only operates on un-instantiated types
"""
    operation_symbol(::AbstractCompositeOperator)
    operation_symbol(::Type{<:AbstractCompositeOperator})

Obtain the model symbol from a given [`AbstractCompositeOperator`](@ref).
"""
operation_symbol(::Type{<:AdditionOperator}) = :(+)
operation_symbol(::Type{<:MultiplicationOperator}) = :(*)
operation_symbol(::Type{<:ConvolutionOperator}) = nothing
operation_symbol(::O) where {O<:AbstractCompositeOperator} = operation_symbol(O)

"""
    CompositeModel{M1,M2,O} <: AbstractSpectralModel
    CompositeModel(left_model, right_model, op::AbstractCompositeOperator)

Type resulting from operations combining any number of [`AbstractSpectralModel`](@ref) via the
model algebra defined from [`AbstractSpectralModelKind`](@ref).

Each operation binary operation in the model algebra is encoded in the parametric types of the
`CompositeModel`, where the operation is given by an [`AbstractCompositeOperator`](@ref).
Composite models adopt the model kind of the `right` model, i.e. `M2`, and obey the model algebra
accordingly.

Composite models very rarely need to be constructed directly, and are instead obtained by regular
model operations.

# Example

```julia
model = XS_PhotoelectricAbsorption() * (XS_PowerLaw() + XS_BlackBody())
typeof(model) <: CompositeModel # true
```
"""
struct CompositeModel{M1,M2,O,T,K} <: AbstractSpectralModel{T,K}
    left::M1
    right::M2
    op::O
    CompositeModel(
        m1::M1,
        m2::M2,
        op::O,
    ) where {M1,M2<:AbstractSpectralModel{T,K},O} where {T,K} = new{M1,M2,O,T,K}(m1, m2, op)
end

modelkind(::Type{<:CompositeModel{M1,M2,O,T,K}}) where {M1,M2,O,T,K} = K()
function implementation(::Type{<:CompositeModel{M1,M2}}) where {M1,M2}
    if (implementation(M1) isa XSPECImplementation) ||
       (implementation(M2) isa XSPECImplementation)
        XSPECImplementation()
    else
        JuliaImplementation()
    end
end
function closurekind(::Type{<:CompositeModel{M1,M2}}) where {M1,M2}
    if (closurekind(M1) isa WithClosures) || (closurekind(M2) isa WithClosures)
        WithClosures()
    else
        WithoutClosures()
    end
end

# invocation wrappers
function invokemodel(e, m::CompositeModel)
    fluxes = construct_objective_cache(m, e)
    invokemodel!(fluxes, e, m)
    view(fluxes, :, 1)
end

function invokemodel!(f, e, model::CompositeModel)
    @assert size(f, 2) == objective_cache_count(model) "Too few flux arrays allocated for this model."
    generated_model_call!(f, e, model, model_parameters_tuple(model))
end

function invokemodel!(f, e, model::CompositeModel, parameters::AbstractArray)
    @assert size(f, 2) == objective_cache_count(model) "Too few flux arrays allocated for this model."
    generated_model_call!(f, e, model, parameters)
end


# algebra grammar
add_models(_, _, ::M1, ::M2) where {M1,M2} =
    throw("Left and right models must be Additive.")
add_models(m1, m2, ::Additive, ::Additive) = CompositeModel(m1, m2, AdditionOperator())
add_models(m1::M1, m2::M2) where {M1,M2} = add_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:+(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel} =
    add_models(m1, m2, modelkind(M1), modelkind(M2))

mult_models(_, _, ::M1, ::M2) where {M1,M2} = throw("Left model must be Multiplicative.")
mult_models(m1, m2, ::Multiplicative, ::AbstractSpectralModelKind) =
    CompositeModel(m1, m2, MultiplicationOperator())
mult_models(m1::M1, m2::M2) where {M1,M2} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:*(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))

conv_models(_, _, ::M1, ::M2) where {M1,M2} =
    throw("Left model must be Convolutional and right model must be Additive.")
conv_models(m1, m2, ::Convolutional, ::Additive) =
    CompositeModel(m1, m2, ConvolutionOperator())
conv_models(m1::M1, m2::M2) where {M1,M2} =
    conv_models(m1, m2, modelkind(M1), modelkind(M2))
(m1::AbstractSpectralModel)(m2::M2) where {M2<:AbstractSpectralModel} = conv_models(m1, m2)

function Base.show(io::IO, @nospecialize(model::CompositeModel))
    expr, infos = _destructure_for_printing(model)
    for (symbol, (m, _)) in zip(keys(infos), infos)
        expr =
            replace(expr, "$(symbol)" => "$(FunctionGeneration.model_base_name(typeof(m)))")
    end
    print(
        io,
        "CompositeModel[",
        Crayons.Crayon(foreground = :cyan),
        "$(expr)",
        Crayons.Crayon(reset = true),
        "]",
    )
end

function _print_param(io, free, name, val, q0, q1, q2, q3, q4)
    print(io, lpad("$name", q0), " ->")
    if val isa FitParam
        info = get_info_tuple(val)
        print(io, lpad(info[1], q1 + 1))
        if free
            print(io, " ± ", rpad(info[2], q2))
            print(io, " ∈ [", lpad(info[3], q3), ", ", rpad(info[4], q4), "]")
            print(
                io,
                Crayons.Crayon(foreground = :green),
                lpad("FREE", 7),
                Crayons.Crayon(reset = true),
            )
        else
            print(
                io,
                Crayons.Crayon(foreground = :cyan),
                lpad("FROZEN", 15 + q1 + q2 + q3 + q4),
                Crayons.Crayon(reset = true),
            )
        end
    end
    println(io)
end

# such an ugly function
function _printinfo(io::IO, model::CompositeModel{M1,M2}) where {M1,M2}
    l_buffer = 5
    expr, infos = _destructure_for_printing(model)
    n_components = length(infos)

    print(io, "CompositeModel with $n_components component models:\n")
    println(
        io,
        Crayons.Crayon(foreground = :cyan),
        " "^l_buffer,
        expr,
        Crayons.Crayon(reset = true),
    )

    info_tuples = reduce(vcat, [get_info_tuple.(modelparameters(i[1])) for i in infos])
    q1, q2, q3, q4 = map(1:4) do i
        maximum(j -> length("$(j[i])"), info_tuples) + 1
    end

    println(io, "Model key and parameters:")
    sym_buffer = 5
    param_name_offset = sym_buffer + maximum(infos) do (_, syms)
        maximum(length(s) for s in syms)
    end
    buff = IOBuffer()
    for (symbol, (m, param_symbols)) in zip(keys(infos), infos)
        M = typeof(m)
        basename = FunctionGeneration.model_base_name(M)
        println(
            buff,
            Crayons.Crayon(foreground = :cyan),
            lpad("$symbol", sym_buffer),
            Crayons.Crayon(reset = true),
            " => ",
            Crayons.Crayon(foreground = :cyan),
            "$basename",
            Crayons.Crayon(reset = true),
        )

        for (val, s::String) in zip(modelparameters(m), param_symbols)
            free = !isfrozen(val)
            _print_param(buff, free, s, val, param_name_offset, q1, q2, q3, q4)
        end
    end
    print(io, String(take!(buff)))
end

function remake_with_number_type(model::CompositeModel)
    error("Todo")
end
# explicitly disable interface for ConstructionBase.jl
ConstructionBase.setproperties(::CompositeModel, ::NamedTuple) =
    throw("Cannot be used with `CompositeModel`.")
ConstructionBase.constructorof(::Type{<:CompositeModel}) =
    throw("Cannot be used with `CompositeModel`.")

function Base.propertynames(model::CompositeModel, private::Bool = false)
    all_parameter_symbols(model)
end

function Base.getproperty(model::CompositeModel, symb::Symbol)
    lookup = all_parameters_to_named_tuple(model)
    lookup[symb]
end

function Base.setproperty!(model::CompositeModel, symb::Symbol, value::FitParam)
    set!(getproperty(model, symb), value)
end

function Base.setproperty!(model::CompositeModel, symb::Symbol, x)
    error(
        "Only `FitParam` may be directly set with another `FitParam`. Use `set_value!` and related API otherwise.",
    )
end

# function ConstructionBase.setproperties(m::CompositeModel, patch::NamedTuple)
# end
