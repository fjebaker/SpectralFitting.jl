export AbstractCompositeOperator,
    AdditionOperator,
    MultiplicationOperator,
    ConvolutionOperator,
    operation_symbol,
    CompositeModel

"""
    abstract type AbstractCompositeOperator

Superype of all composition operators. Used to implement the model algebra of
[`AbstractSpectralModelKind`](@ref) through a trait system.

The Julia symbol corresponding to a given `AbstractCompositeOperator` may be
obtained through [`operation_symbol`](@ref).
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

"""
    operation_symbol(::AbstractCompositeOperator)
    operation_symbol(::Type{<:AbstractCompositeOperator})

Obtain the model symbol from a given [`AbstractCompositeOperator`](@ref).
"""
operation_symbol(::Type{<:AdditionOperator}) = :(+)
operation_symbol(::Type{<:MultiplicationOperator}) = :(*)
operation_symbol(::O) where {O<:AbstractCompositeOperator} = operation_symbol(O)

combine_components!(::Type{<:AdditionOperator}, left, right) = @. right = left + right
combine_components!(::Type{<:MultiplicationOperator}, left, right) = @. right = left * right

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

"""
    CompositeModel{T,K,Op} <: AbstractSpectralModel{T,K}
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
model = PhotoelectricAbsorption() * (PowerLaw() + BlackBody())
typeof(model) <: CompositeModel # true
```
"""
struct CompositeModel{T,K,Operator,M1<:AbstractSpectralModel,M2<:AbstractSpectralModel} <:
       AbstractSpectralModel{T,K}
    left::M1
    right::M2
end

function CompositeModel(
    left::AbstractSpectralModel{T},
    right::AbstractSpectralModel{T,K},
    ::Operator,
) where {T,K,Operator}
    CompositeModel{T,K,Operator,typeof(left),typeof(right)}(copy(left), copy(right))
end

function Base.copy(m::CompositeModel{T,K,Op}) where {T,K,Op}
    CompositeModel(getfield(m, :left), getfield(m, :right), Op())
end

function implementation(::Type{<:CompositeModel{T,K,Op,M1,M2}}) where {T,K,Op,M1,M2}
    if (implementation(M1) isa XSPECImplementation) ||
       (implementation(M2) isa XSPECImplementation)
        XSPECImplementation()
    else
        JuliaImplementation()
    end
end

objective_cache_count(model::AbstractSpectralModel) = 1
objective_cache_count(model::CompositeModel) =
    objective_cache_count(getfield(model, :left)) +
    objective_cache_count(getfield(model, :right))

parameter_vector(model::CompositeModel) = vcat(
    parameter_vector(getfield(model, :left)),
    parameter_vector(getfield(model, :right)),
)
parameter_count(model::CompositeModel) =
    parameter_count(getfield(model, :left)) + parameter_count(getfield(model, :right))

function remake_with_number_type(model::CompositeModel{T,K,Op}) where {T,K,Op}
    CompositeModel(
        remake_with_number_type(getfield(model, :left)),
        remake_with_number_type(getfield(model, :right)),
        Op(),
    )
end

function remake_with_parameters(model::CompositeModel{T,K,Op}, params::Tuple) where {T,K,Op}
    P_left = parameter_count(getfield(model, :left))
    CompositeModel(
        remake_with_parameters(getfield(model, :left), params[1:P_left]),
        remake_with_parameters(getfield(model, :right), params[(P_left+1):end]),
        Op(),
    )
end

# invocation wrappers
function invokemodel(e, m::CompositeModel)
    fluxes = allocate_model_output(m, e)
    invokemodel!(fluxes, e, m)
    view(fluxes, :, 1)
end

function invokemodel!(
    outputs,
    domain,
    model::CompositeModel{T,K,Op,M1,M2},
) where {T,K,Op,M1,M2}
    @assert size(outputs, 2) >= objective_cache_count(model) "Too few flux arrays allocated for this model ($(size(outputs)))."
    RetType = typeof(@views(outputs[:, 1]))

    @views begin
        _out_right = if M2 <: CompositeModel
            outputs[:, 1:end]
        else
            outputs[:, 1:1]
        end

        invokemodel!(_out_right, domain, getfield(model, :right))
        if Op <: ConvolutionOperator
            # Use the same principle output as the right hand side
            _out_left = if M1 <: CompositeModel
                outputs[:, 1:end]
            else
                outputs[:, 1:1]
            end
            invokemodel!(_out_left, domain, getfield(model, :left))
        else
            _out_left = if M1 <: CompositeModel
                outputs[:, 2:end]
            else
                outputs[:, 2:2]
            end
            invokemodel!(_out_left, domain, getfield(model, :left))
            combine_components!(Op, outputs[:, 2], outputs[:, 1])
        end
    end
end

# TODO
# function invokemodel!(
#     outputs,
#     domain,
#     model::CompositeModel{T,K,Op},
#     parameters::AbstractVector,
# ) where {T,K,Op}
#     @assert size(outputs, 2) >= objective_cache_count(model) "Too few flux arrays allocated for this model."
#     P_left = parameter_count(getfield(model, :left))

#     @views begin
#         invokemodel!(
#             outputs[:, 1:end],
#             domain,
#             getfield(model, :right),
#             parameters[(P_left+1):end],
#         )

#         if Op <: ConvolutionOperator
#             invokemodel!(
#                 outputs[:, 1:end],
#                 domain,
#                 getfield(model, :left),
#                 parameters[1:P_left],
#             )
#         else
#             invokemodel!(
#                 outputs[:, 2:end],
#                 domain,
#                 getfield(model, :left),
#                 parameters[1:P_left],
#             )
#             combine_components!(Op, outputs[:, 2], outputs[:, 1])
#         end
#     end
# end

# printing
struct ModelKindCounts{N}
    additive::Int
    multiplicative::Int
    convolutional::Int
    syms::NTuple{N,Symbol}
end

function destructure(::Type{<:AbstractSpectralModel{T,K}}, d::ModelKindCounts) where {T,K}
    if K == Additive
        s = Symbol("a$(d.additive + 1)")
        s,
        ModelKindCounts(d.additive + 1, d.multiplicative, d.convolutional, (d.syms..., s))
    elseif K == Multiplicative
        s = Symbol("m$(d.multiplicative + 1)")
        s,
        ModelKindCounts(d.additive, d.multiplicative + 1, d.convolutional, (d.syms..., s))
    elseif K == Convolutional
        s = Symbol("c$(d.convolutional + 1)")
        s,
        ModelKindCounts(d.additive, d.multiplicative, d.convolutional + 1, (d.syms..., s))
    else
        error("Unknown model kind: $K")
    end
end

function destructure(
    ::Type{<:CompositeModel{T,K,Op,Left,Right}},
    d::ModelKindCounts,
) where {T,K,Op,Left,Right}
    exp_left, d_left = destructure(Left, d)
    exp_right, d_right = destructure(Right, d_left)
    if Op <: ConvolutionOperator
        :($(exp_left)($exp_right)), d_right
    else
        op = SpectralFitting.operation_symbol(Op)
        :($(op)($exp_left, $exp_right)), d_right
    end
end

destructure(M::Type{<:CompositeModel}) = destructure(M, ModelKindCounts(0, 0, 0, ()))
flatten_models(m::AbstractSpectralModel) = (m,)

function flatten_models(m::CompositeModel)
    right = flatten_models(getfield(m, :right))
    left = flatten_models(getfield(m, :left))
    (left..., right...)
end

struct DestructuredCompositeModel{M}
    expr::String
    models::Vector{Pair{Symbol,M}}
end

function destructure(@nospecialize(m::CompositeModel))
    expr, info = destructure(typeof(m))
    all_models = flatten_models(m)
    DestructuredCompositeModel(
        "$(expr)",
        Pair{Symbol,AbstractSpectralModel}[k => v for (k, v) in zip(info.syms, all_models)],
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
        end

        if free
            printstyled(io, lpad("FREE", 7), color = :green)
        else
            print(io, "  ")
            printstyled(io, lpad("FROZEN", 15 + q1 + q2 + q3 + q4), color = :cyan)
        end
    else
        print(io, val)
    end
    println(io)
end


function _printinfo(io::IO, @nospecialize(model::CompositeModel); bindings = nothing)
    expr_buffer = 3
    sym_buffer = 6

    dis = destructure(model)

    println(io, "CompositeModel with $(length(dis.models)) model components:")
    print(io, " "^expr_buffer)
    printstyled(io, dis.expr, color = :cyan)
    println(io)
    println(io, "Model key and parameters:")

    param_index = 1
    for (sym, m) in dis.models
        basename = Base.typename(typeof(m)).name

        printstyled(io, lpad("$sym", sym_buffer), color = :cyan)
        print(io, " => ")

        buff = IOBuffer()
        _printinfo(IOContext(buff, io), m; bindings = if isnothing(bindings)
            nothing
        else
            get(bindings, sym, nothing)
        end)
        s = String(take!(buff))
        println(io, strip(indent(s, 4)))
    end
end

# explicitly disable interface for ConstructionBase.jl
ConstructionBase.setproperties(::CompositeModel, ::NamedTuple) =
    throw("Cannot be used with `CompositeModel`.")
ConstructionBase.constructorof(::Type{<:CompositeModel}) =
    throw("Cannot be used with `CompositeModel`.")


function Base.propertynames(model::CompositeModel)
    _, info = destructure(typeof(model))
    info.syms
end

function Base.getproperty(model::CompositeModel, symb::Symbol)
    model_map = destructure(model)
    ind = findfirst(i -> i[1] == symb, model_map.models)
    if isnothing(ind)
        error("$(_model_name(model)) has no component $symb")
    else
        last(model_map.models[ind])
    end
end
