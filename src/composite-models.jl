export AbstractCompositeOperator,
    AdditionOperator,
    MultiplicationOperator,
    ConvolutionOperator,
    operation_symbol,
    CompositeModel,
    get_free_model_params,
    get_frozen_model_params

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
struct CompositeModel{M1,M2,O} <: AbstractSpectralModel
    left::M1
    right::M2
    op::O
end

modelkind(::Type{<:CompositeModel{M1,M2}}) where {M1,M2} = modelkind(M2)
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

# parameter access
# get_param_symbols(M::Type{<:AbstractSpectralModel}) = fieldnames(M)
function get_param_types(model::Type{<:SpectralFitting.CompositeModel})
    __generated_get_param_types(model)
end

# utilities
add_param_types!(_, ::Type{Nothing}) = nothing
add_param_types!(types, M::Type{<:AbstractSpectralModel}) =
    push!(types, get_param_types(M)...)

# invokation wrappers
function invokemodel!(f, e, model::CompositeModel)
    params = get_params_value(model)
    generated_model_call!(f, e, model, params)
end

function invokemodel(e, m::CompositeModel)
    fluxes = make_fluxes(e, flux_count(m))
    invokemodel!(fluxes, e, m)
    first(fluxes)
end
function invokemodel(e, m::CompositeModel, free_params)
    if eltype(free_params) <: Number
        # for compatability with AD
        fluxes = make_fluxes(e, flux_count(m), eltype(free_params))
        invokemodel!(fluxes, e, m, free_params)
    else
        p0 = get_value.(free_params)
        fluxes = make_fluxes(e, flux_count(m), eltype(p0))
        invokemodel!(fluxes, e, m, p0)
    end
    first(fluxes)
end
function invokemodel!(f, e, m::CompositeModel, free_params)
    frozen_params = get_value.(get_frozen_model_params(m))
    invokemodel!(f, e, m, free_params, frozen_params)
end
function invokemodel!(f, e, model::CompositeModel, free_params, frozen_params)
    generated_model_call!(f, e, model, free_params, frozen_params)
end

# algebra grammar
add_models(_, _, ::M1, ::M2) where {M1,M2} =
    error("Left and right models must be Additive.")
add_models(m1, m2, ::Additive, ::Additive) =
    CompositeModel(m1, m2, AdditionOperator())
add_models(m1::M1, m2::M2) where {M1,M2} = add_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:+(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel} =
    add_models(m1, m2, modelkind(M1), modelkind(M2))

mult_models(_, _, ::M1, ::M2) where {M1,M2} = error("Left model must be Multiplicative.")
mult_models(m1, m2, ::Multiplicative, ::AbstractSpectralModelKind) =
    CompositeModel(m1, m2, MultiplicationOperator())
mult_models(m1::M1, m2::M2) where {M1,M2} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:*(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))

conv_models(_, _, ::M1, ::M2) where {M1,M2} =
    error("Left model must be Convolutional and right model must be Additive.")
conv_models(m1, m2, ::Convolutional, ::Additive) =
    CompositeModel(m1, m2, ConvolutionOperator())
conv_models(m1::M1, m2::M2) where {M1,M2} =
    conv_models(m1, m2, modelkind(M1), modelkind(M2))
(m1::AbstractSpectralModel)(m2::M2) where {M2<:AbstractSpectralModel} = conv_models(m1, m2)

# runtime param access
function _add_params_to_index!(params, adder!, m)
    if !isnothing(m)
        foreach(get_params(m)) do p
            adder!(params, p)
        end
    end
end

__get_model_parameters!(params, adder!, model::AbstractSpectralModel) =
    _add_params_to_index!(params, adder!, model)

function __get_model_parameters!(params, adder!, model::CompositeModel)
    recursive_model_parse(model) do (left, right, _)
        _add_params_to_index!(params, adder!, right)
        _add_params_to_index!(params, adder!, left)
        nothing
    end
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

function get_params_value(model::CompositeModel)
    T = generated_get_model_number_type(model)
    params = T[]
    __get_model_parameters!(params, (ps, p) -> push!(ps, get_value(p)), model)
    params
end

function get_params(model::CompositeModel)
    params = AbstractFitParameter[]
    __get_model_parameters!(params, push!, model)
    params
end

function _add_symbols_and_params_to_index!(params, m)
    if !isnothing(m)
        foreach(get_param_symbol_pairs(m)) do (s, p)
            symb = _make_unique_readable_symbol(s, first.(params))
            push!(params, symb => p)
        end
    end
end

function get_param_symbol_pairs(model::CompositeModel)
    params = Pair{Symbol,AbstractFitParameter}[]
    recursive_model_parse(model) do (left, right, _)
        _add_symbols_and_params_to_index!(params, right)
        _add_symbols_and_params_to_index!(params, left)
        nothing
    end
    params
end

function get_param(model::CompositeModel, s::Symbol)
    ps = get_param_symbol_pairs(model)
    i = findfirst(i -> first(i) == s, ps)
    if !isnothing(i)
        last(ps[i])
    else
        error("Model has no symbol $s.")
    end
end

# printing
_expression_string(model::M) where {M<:AbstractSpectralModel} = modelinfo(M)
_expression_string(left, right, ::Multiplicative) = "$left * $right"
_expression_string(left, right, ::Convolutional) = "$left($right)"
_expression_string(left, right, ::Additive) = "($left + $right)"

function _expression_string(cm::CompositeModel{M1,M2}) where {M1,M2}
    left = _expression_string(cm.left)
    right = _expression_string(cm.right)
    _expression_string(left, right, modelkind(M1))
end

function Base.show(io::IO, ::CompositeModel)
    print(io, "CompositeModel")
end

# such an ugly function
function _printinfo(io::IO, model::CompositeModel{M1,M2}) where {M1,M2}
    l_buffer = 5
    model_types = generated_get_model_types(model)
    n_components = length(model_types)
    expr, infos = _readable_expression_info(model)

    print(io, "CompositeModel with $n_components component models:\n")
    println(
        io,
        Crayons.Crayon(foreground = :cyan),
        " "^l_buffer,
        "$expr",
        Crayons.Crayon(reset = true),
    )

    params = Pair{Symbol,AbstractFitParameter}[]
    model_list = map(infos.models) do (symbol, m)
        _ps = map(get_param_symbol_pairs(m)) do (ps, p)
            ups = _make_unique_readable_symbol(ps, first.(params))
            push!(params, ups => p)
            String(ups)
        end
        p_string = join(_ps, ", ")
        symbol => ("$(model_base_name(typeof(m)))", p_string)
    end

    pad1 = maximum(i -> length(String(i[1])), params) + l_buffer
    pad2 = maximum(i -> length(i[2][1]), model_list) + 1

    println(io, lpad("Symbols:", pad1), "    ", rpad("Models:", pad2), " Params:")
    for (s, (m, mp)) in model_list
        println(io, "$(lpad(s, pad1)) => $(rpad(m, pad2)) ($mp) ")
    end

    println(io, lpad("Symbols:", pad1), "     Values:")

    info_tuple = get_info_tuple.(last.(params))
    q1, q2, q3, q4 = map(1:4) do i
        maximum(j -> length("$(j[i])"), info_tuple) + 1
    end
    for ((s, p), info) in zip(params, info_tuple)
        print(io, "$(lpad(s, pad1)) => ")
        print(io, lpad(info[1], q1))
        if !is_frozen(p)
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
        println(io)
    end
end

_unique_model_symbol(::M, infos) where {M<:AbstractSpectralModel} =
    _unique_model_symbol(modelkind(M), infos)
_unique_model_symbol(::Additive, infos) = Symbol('a', infos.a[] += 1)
_unique_model_symbol(::Multiplicative, infos) = Symbol('m', infos.m[] += 1)
_unique_model_symbol(::Convolutional, infos) = Symbol('c', infos.c[] += 1)

function _readable_expression_info(model::CompositeModel)
    models = Pair{Symbol,AbstractSpectralModel}[]
    infos = (models = models, a = Ref(0), m = Ref(0), c = Ref(0))
    expr = recursive_model_parse(model) do (left, right, op)
        r_symb = if right isa AbstractSpectralModel
            rs = _unique_model_symbol(right, infos)
            push!(models, rs => right)
            rs
        else
            right
        end
        l_symb = if left isa AbstractSpectralModel
            ls = _unique_model_symbol(left, infos)
            push!(models, ls => left)
            ls
        else
            left
        end
        if (op == :*)
            :($(l_symb) * $(r_symb))
        elseif (op == :+)
            :($(l_symb) + $(r_symb))
        else
            :($(l_symb)($r_symb))
        end
    end
    expr, infos
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    cm::CompositeModel{M1,M2},
) where {M1,M2}
    _printinfo(io, cm)
end
