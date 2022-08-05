# model parameter introspection

parameter_symbols_and_types(M::Type{<:AbstractSpectralModel}) = (fieldnames(M), M.types)

push_symbol_and_type!(_, _, ::Type{Nothing}) = nothing
function push_symbol_and_type!(names, types, ::Type{M}) where {M<:AbstractSpectralModel}
    for (s, t) in zip(parameter_symbols_and_types(M)...)
        push!(names, s)
        push!(types, t)
    end
end

function model_parameter_symbol_and_type(
    model::Type{<:SpectralFitting.CompositeSpectralModel},
)
    names = Symbol[]
    types = Type[]
    SpectralFitting.recursive_model_parse(model) do (left, right, _)
        push_symbol_and_type!(names, types, right)
        push_symbol_and_type!(names, types, left)
    end
    names, types
end

function model_parameter_symbol_and_type(model::Type{<:AbstractSpectralModel})
    names = Symbol[]
    types = Type[]
    push_symbol_and_type!(names, types, model)
    names, types
end

function model_parameter_count(M::Type{<:AbstractSpectralModel})::Int
    _, types = model_parameter_symbol_and_type(M)
    length(types)
end

function model_frozen_parameter_count(M::Type{<:AbstractSpectralModel})::Int
    _, types = model_parameter_symbol_and_type(M)
    count(is_frozen, types)
end

function model_free_parameter_count(M::Type{<:AbstractSpectralModel})::Int
    _, types = model_parameter_symbol_and_type(M)
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
    parameters,
    flux,
    ::Type{M},
) where {M<:AbstractSpectralModel}
    params = Symbol[]

    for p in fieldnames(M)
        i = Base.gensym(p)
        push!(parameters, i)
        push!(params, i)
    end

    base_name = Base.typename(M).name
    :(invokemodel!($flux, energy, $(base_name){Float64}, $(params...)))
end

function add_flux_invokation!(flux_count, parameters, M, K::Convolutional)
    flux = Symbol(:flux, flux_count[])
    assemble_invoke_statment!(parameters, flux, M)
end

function add_flux_invokation!(flux_count, parameters, M, K)
    i = (flux_count[] += 1)
    flux = Symbol(:flux, i)
    assemble_invoke_statment!(parameters, flux, M)
end

function add_flux_invokation!(
    flux_count,
    parameters,
    ::Type{M},
) where {M<:AbstractSpectralModel}
    add_flux_invokation!(flux_count, parameters, M, modelkind(M))
end

push_resolve_flux_combine!(_, _, ::Type{ConvolutionOperator}) = nothing
function push_resolve_flux_combine!(
    statements,
    flux_count,
    ::Type{O},
) where {O<:AbstractCompositeOperator}
    push_resolve_flux_combine!(statements, flux_count, operation_symbol(O))
end

function push_resolve_flux_combine!(statements, flux_count, op::Symbol)
    i = (flux_count[] -= 1)
    flux_right = Symbol(:flux, i + 1)
    flux_left = Symbol(:flux, i)
    expr = Expr(:call, op, flux_left, flux_right)
    push!(statements, :(@.($flux_left = $expr)))
end

push_model!(_, _, _, ::Type{Nothing}) = nothing
function push_model!(
    statements,
    parameters,
    flux_count,
    model::Type{<:AbstractSpectralModel},
)
    push!(statements, add_flux_invokation!(flux_count, parameters, model))
end

__adjust_max!(flux_max, flux_counts) = flux_max[] = (flux_counts[] > flux_max[] ? flux_counts[] : flux_max[])
function assemble_expression(model::Type{<:CompositeSpectralModel})
    flux_count = Ref(0)
    flux_count_max = Ref(0)
    statements = Expr[]
    parameters = Symbol[]
    recursive_model_parse(model) do (left, right, op)
        push_model!(statements, parameters, flux_count, right)
        __adjust_max!(flux_count_max, flux_count)
        push_model!(statements, parameters, flux_count, left)
        __adjust_max!(flux_count_max, flux_count)
        push_resolve_flux_combine!(statements, flux_count, op)
        __adjust_max!(flux_count_max, flux_count)
    end
    statements, parameters, flux_count_max[]
end

"""
    assemble_expression(model::Type{<:AbstractSpectralModel})

Specialisation for when only a single model invoked.
"""
function assemble_expression(model::Type{<:AbstractSpectralModel})
    parameters = Symbol[]
    s = assemble_invoke_statment!(parameters, :flux1, model)
    (s,), parameters
end

function maximum_flux_count(model)
    _, _, fc = assemble_expression(model)
    fc
end

@generated function generated_maximum_flux_count(model)
    N = maximum_flux_count(model)
    :($N)
end

@generated function generated_model_call!(fluxes, energy, model, params)
    expr, parameters, _ = assemble_expression(model)
    flux_unpack = [Symbol(:flux, i) for i in 1:maximum_flux_count(model)]
    p_assign = [:($s = params[$i]) for (i, s) in enumerate(parameters)]
    quote
        @fastmath begin
            @inbounds let ($(flux_unpack...),) = fluxes
                $(p_assign...)
                $(expr...)
                return flux1
            end
        end
    end
end

@generated function generated_model_call!(fluxes, energy, model, free_params, frozen_params)
    expr, parameters, _ = assemble_expression(model)
    _, p_types = model_parameter_symbol_and_type(model)

    flux_unpack = [Symbol(:flux, i) for i in 1:maximum_flux_count(model)]

    # unpack free and frozen seperately
    i_frozen = 0
    i_free = 0
    p_assign = map(enumerate(parameters)) do (i, s)
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
                $(expr...)
                return flux1
            end
        end
    end
end
