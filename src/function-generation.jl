# assemble_invoke_statment!(params, flux, M)
# parameter_symbols_and_types(M)

mutable struct GenerationAggragate
    statements::Vector{Expr}
    parameters::Vector{Symbol}
    closure_parameters::Vector{Symbol}
    models::Vector{Type}
    flux_count::Int
    maximum_flux_count::Int
    GenerationAggragate() = new(Expr[], Symbol[], Symbol[], Type[], 0, 0)
end

push_closure_param!(g::GenerationAggragate, s::Symbol) = push!(g.closure_parameters, s)
push_param!(g::GenerationAggragate, s::Symbol) = push!(g.parameters, s)
push_statement!(g::GenerationAggragate, s::Expr) = push!(g.statements, s)
push_model!(g::GenerationAggragate, t::Type) = push!(g.models, t)

function new_parameter!(ga::GenerationAggragate, p::Symbol)
    param = Base.gensym(p)
    push_param!(ga, param)
    param
end

get_flux_symbol(i::Int) = Symbol(:flux, i)
function get_flux_symbol!(g::GenerationAggragate)
    i = inc_flux!(g)
    Symbol(:flux, i)
end

function set_flux!(g::GenerationAggragate, f)
    g.flux_count = f
    if g.flux_count > g.maximum_flux_count
        g.maximum_flux_count = g.flux_count
    end
    g.flux_count
end

inc_flux!(g::GenerationAggragate) = set_flux!(g, g.flux_count + 1)
dec_flux!(g::GenerationAggragate) = set_flux!(g, g.flux_count - 1)


# model parameter introspection

model_parameter_types(M::Type{<:AbstractSpectralModel}) = M.types

push_parameter_types!(_, ::Type{Nothing}) = nothing
function push_parameter_types!(types, M::Type{<:AbstractSpectralModel})
    for t in model_parameter_types(M)
        push!(types, t)
    end
end

function model_parameter_types(
    model::Type{<:SpectralFitting.CompositeSpectralModel},
)
    types = Type[]
    SpectralFitting.recursive_model_parse(model) do (left, right, _)
        push_parameter_types!(types, right)
        push_parameter_types!(types, left)
    end
    types
end

function model_parameter_count(M::Type{<:AbstractSpectralModel})::Int
    types = model_parameter_types(M)
    length(types)
end

function model_frozen_parameter_count(M::Type{<:AbstractSpectralModel})::Int
    types = model_parameter_types(M)
    count(is_frozen, types)
end

function model_free_parameter_count(M::Type{<:AbstractSpectralModel})::Int
    types = model_parameter_types(M)
    count(!is_frozen, types)
end

@generated function generated_model_parameter_count(model)
    N = model_parameter_count(model)
    :($N)
end

@generated function generated_frozen_model_parameter_count(model)
    N = model_frozen_parameter_count(model)
    :($N)
end

@generated function generated_free_model_parameter_count(model)
    N = model_free_parameter_count(model)
    :($N)
end

# model function assembly

function assemble_invoke_statment!(
    ga::GenerationAggragate,
    flux,
    M::Type{<:AbstractSpectralModel},
)
    params = map(fieldnames(M)) do p
        new_parameter!(ga, p)
    end

    base_name = Base.typename(M).name
    :(invokemodel!($flux, energy, $(base_name), $(params...)))
end

function add_flux_invokation!(ga::GenerationAggragate, M, ::Convolutional)
    flux = get_flux_symbol(ga.flux_count)
    s = assemble_invoke_statment!(ga, flux, M)
    push_statement!(ga, s)
end

function add_flux_invokation!(ga::GenerationAggragate, M, ::AbstractSpectralModelKind)
    inc_flux!(ga)
    flux = get_flux_symbol(ga.flux_count)
    s = assemble_invoke_statment!(ga, flux, M)
    push_statement!(ga, s)
end

add_flux_invokation!(ga::GenerationAggragate, M::Type{<:AbstractSpectralModel}) =
    add_flux_invokation!(ga, M, modelkind(M))

assemble_and_add_resolve_flux_combine!(::GenerationAggragate, ::Type{ConvolutionOperator}) = nothing
assemble_and_add_resolve_flux_combine!( ga::GenerationAggragate, O::Type{<:AbstractCompositeOperator}) =
    assemble_and_add_resolve_flux_combine!(ga, operation_symbol(O))

function assemble_and_add_resolve_flux_combine!(ga::GenerationAggragate, op::Symbol)
    dec_flux!(ga)
    fr = get_flux_symbol(ga.flux_count + 1)
    fl = get_flux_symbol(ga.flux_count)
    expr = Expr(:call, op, fl, fr)
    push_statement!(ga, :(@.($fl = $expr)))
end

assemble_and_add_model!(::GenerationAggragate, ::Type{Nothing}) = nothing
assemble_and_add_model!(ga::GenerationAggragate, model::Type{<:AbstractSpectralModel}) =
    add_flux_invokation!(ga, model)

__adjust_max!(flux_max, flux_counts) = flux_max[] = (flux_counts[] > flux_max[] ? flux_counts[] : flux_max[])
function assemble_expression(model::Type{<:CompositeSpectralModel})
    ga = GenerationAggragate()
    recursive_model_parse(model) do (left, right, op)
        assemble_and_add_model!(ga, right)
        assemble_and_add_model!(ga, left)
        assemble_and_add_resolve_flux_combine!(ga, op)
    end
    return ga
end

"""
    assemble_expression(model::Type{<:AbstractSpectralModel})

Specialisation for when only a single model invoked.
"""
function assemble_expression(model::Type{<:AbstractSpectralModel})
    ga = GenerationAggragate()
    flux = get_flux_symbol(inc_flux!(ga))
    s = assemble_invoke_statment!(ga, flux, model)
    push_statement!(ga, s)
    ga
end

function maximum_flux_count(model)
    ga = assemble_expression(model)
    ga.maximum_flux_count
end

@generated function generated_maximum_flux_count(model)
    N = maximum_flux_count(model)
    :($N)
end

@generated function generated_model_call!(fluxes, energy, model, params)
    ga = assemble_expression(model)
    flux_unpack = [Symbol(:flux, i) for i in 1:ga.maximum_flux_count]
    p_assign = [:($s = params[$i]) for (i, s) in enumerate(ga.parameters)]
    quote
        @fastmath begin
            @inbounds let ($(flux_unpack...),) = fluxes
                $(p_assign...)
                $(ga.statements...)
                return flux1
            end
        end
    end
end

@generated function generated_model_call!(fluxes, energy, model, free_params, frozen_params)
    ga = assemble_expression(model)
    p_types = model_parameter_types(model)

    flux_unpack = [Symbol(:flux, i) for i in 1:ga.maximum_flux_count]

    # unpack free and frozen seperately
    i_frozen = 0
    i_free = 0
    p_assign = map(enumerate(ga.parameters)) do (i, s)
        if is_frozen(p_types[i])
            :($s = frozen_params[$(i_frozen += 1)])
        else
            :($s = free_params[$(i_free += 1)])
        end
    end

    quote
        @fastmath begin
            @inbounds let ($(flux_unpack...),) = fluxes
                $(p_assign...)
                $(ga.statements...)
                return flux1
            end
        end
    end
end
