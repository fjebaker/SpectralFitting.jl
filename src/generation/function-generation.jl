# add_invoke_statment!(params, flux, M)
# parameter_symbols_and_types(M)

module FunctionGeneration

import SpectralFitting
import SpectralFitting:
    AbstractSpectralModel,
    AbstractSpectralModelKind,
    Convolutional,
    CompositeModel,
    operation_symbol,
    modelkind,
    has_closure_params

include("parsing-utilities.jl")

mutable struct GenerationAggregate{NumType}
    statements::Vector{Expr}
    infos::Vector{ModelInfo}
    closure_params::Vector{Symbol}
    models::Vector{Type}
    flux_count::Int
    maximum_flux_count::Int
end
GenerationAggregate(T, infos) =
    GenerationAggregate{T}(Expr[], [infos], Symbol[], Type[], 0, 0)
GenerationAggregate(T, infos::AbstractVector) =
    GenerationAggregate{T}(Expr[], infos, Symbol[], Type[], 0, 0)

push_closure_param!(g::GenerationAggregate, s::Symbol) = push!(g.closure_params, s)
push_info!(g::GenerationAggregate, s::ModelInfo) = push!(g.infos, s)
push_statement!(g::GenerationAggregate, s::Expr) = push!(g.statements, s)
push_model!(g::GenerationAggregate, t::Type) = push!(g.models, t)

function new_closure_param!(ga::GenerationAggregate, p::Symbol)
    param = Base.gensym(p)
    push_closure_param!(ga, param)
    param
end

get_flux_symbol(i::Int) = Symbol(:flux, i)
function get_flux_symbol!(g::GenerationAggregate)
    i = inc_flux!(g)
    Symbol(:flux, i)
end

function set_flux!(g::GenerationAggregate, f)
    g.flux_count = f
    if g.flux_count > g.maximum_flux_count
        g.maximum_flux_count = g.flux_count
    end
    g.flux_count
end

inc_flux!(g::GenerationAggregate) = set_flux!(g, g.flux_count + 1)
dec_flux!(g::GenerationAggregate) = set_flux!(g, g.flux_count - 1)

# model invokation generation
function add_invoke_statment!(
    ga::GenerationAggregate{NumType},
    flux,
    M::Type{<:AbstractSpectralModel},
) where {NumType}
    index = length(ga.models) + 1
    # infos are in the correct order if indexes this way
    info = ga.infos[index]

    # aggregate closure parameters
    closure_params = if has_closure_params(M)
        map(closure_parameter_symbols(M)) do p
            new_closure_param!(ga, p)
        end
    else
        ()
    end

    # get the parameter type
    T = NumType <: SpectralFitting.FitParam ? NumType.parameters[1] : NumType
    model_constructor =
        :($(M.name.wrapper){$(M.parameters[1:end-2]...),$(T),$(M.parameters[end])})

    # assemble the invocation statement
    s = :(invokemodel!(
        $flux,
        energy,
        $(model_constructor)($(closure_params...), $(info.generated_symbols...)),
    ))
    push_model!(ga, M)
    push_statement!(ga, s)
end

# don't increment flux for Convolutional models
function add_invokation!(ga::GenerationAggregate, M, ::Convolutional)
    flux = get_flux_symbol(ga.flux_count)
    add_invoke_statment!(ga, flux, M)
end
# increment for everything else
function add_invokation!(ga::GenerationAggregate, M, ::AbstractSpectralModelKind)
    inc_flux!(ga)
    flux = get_flux_symbol(ga.flux_count)
    add_invoke_statment!(ga, flux, M)
end

function add_flux_resolution!(ga::GenerationAggregate, op::Symbol)
    fr = get_flux_symbol(ga.flux_count)
    dec_flux!(ga)
    fl = get_flux_symbol(ga.flux_count)
    expr = Expr(:call, op, fl, fr)
    push_statement!(ga, :(@.($fl = $expr)))
end

function assemble_closures(ga::GenerationAggregate, model)
    info_with_closures = filter(ga.infos) do info
        has_closure_params(info.type)
    end
    i = 0
    closures = map(info_with_closures) do info
        cl = map(closure_parameter_symbols(info.type)) do s
            param = ga.closure_params[i+=1]
            path = :(getproperty($(info.lens), $(Meta.quot(s))))
            :($param = $(path))
        end
        [cl...]
    end
    if isempty(closures)
        ()
    else
        reduce(vcat, closures)
    end
end

function generated_maximum_flux_count(model::Type{<:AbstractSpectralModel{T}}) where {T}
    ga = assemble_aggregate_info(model, T)
    :($(ga.maximum_flux_count))
end

function assemble_parameter_assignment(ga::GenerationAggregate, model)
    # unpack free and frozen seperately
    i_frozen = 0
    i_free = 0
    all = map(ga.infos) do info
        assignments = map(zip(info.symbols, info.generated_symbols)) do ((p, s))
            if (p in info.frozen)
                :($(s) = frozen_params[$(i_frozen += 1)])
            else
                :($(s) = free_params[$(i_free += 1)])
            end
        end
        assignments
    end
    reduce(vcat, all)
end

function generate_call(flux_unpack, closures, p_assignments, statements)
    quote
        @fastmath begin
            @inbounds let ($(flux_unpack...),) = fluxes
                $(closures...)
                $(p_assignments...)
                $(statements...)
                return flux1
            end
        end
    end
end

function generated_model_call!(fluxes, energy, model, free_params, frozen_params)
    # propagate information about free parameters to allow for AD
    ga = assemble_aggregate_info(model, eltype(free_params))
    flux_unpack = [Symbol(:flux, i) for i = 1:ga.maximum_flux_count]
    p_assign = assemble_parameter_assignment(ga, model)
    closures = assemble_closures(ga, model)
    generate_call(flux_unpack, closures, p_assign, ga.statements)
end
function generated_model_call!(fluxes, energy, model, params)
    # propagate information about free parameters to allow for AD
    ga = assemble_aggregate_info(model, eltype(params))
    flux_unpack = [Symbol(:flux, i) for i = 1:ga.maximum_flux_count]
    i = 0
    p_assign = reduce(
        vcat,
        [
            [:($(s) = params[$(i += 1)]) for s in info.generated_symbols] for
            info in ga.infos
        ],
    )
    closures = assemble_closures(ga, model)
    generate_call(flux_unpack, closures, p_assign, ga.statements)
end


function assemble_aggregate_info(model::Type{<:CompositeModel}, NumType)
    ga = FunctionGeneration.GenerationAggregate(NumType, getinfo(model))
    FunctionGeneration.recursive_model_parse(model) do (left, right, op_type)
        # get operation symbol
        op = operation_symbol(op_type)
        if (right !== Nothing)
            FunctionGeneration.add_invokation!(ga, right, modelkind(right))
        end
        if (left !== Nothing)
            FunctionGeneration.add_invokation!(ga, left, modelkind(left))
        end
        if (!isnothing(op))
            FunctionGeneration.add_flux_resolution!(ga, op)
        end
        Nothing
    end
    ga
end
function assemble_aggregate_info(model::Type{<:AbstractSpectralModel}, NumType)
    ga = FunctionGeneration.GenerationAggregate(NumType, getinfo(model))
    flux = FunctionGeneration.get_flux_symbol(FunctionGeneration.inc_flux!(ga))
    FunctionGeneration.add_invoke_statment!(ga, flux, model)
    ga
end

model_T(::Type{<:AbstractSpectralModel{T}}) where {T} = T

function rebuild_composite_model(model)
    expr = recursive_model_parse(model) do (left, right, op)
        op_symbol = operation_symbol(op)
        Expr(:call, op_symbol, :($(left)()), :($(right)()))
    end
    quote
        $(expr)
    end
end

function modelparameterstuple(model::Type{<:AbstractSpectralModel})
    info = getinfo(model)
    _parameter_lenses(info, info.symbols)
end
function modelparameterstuple(model::Type{<:CompositeModel})
    infos = getinfo(model)
    reduce(vcat, map(i -> _parameter_lenses(i, i.symbols), infos))
end

function freeparameterstuple(model::Type{<:AbstractSpectralModel})
    info = getinfo(model)
    _parameter_lenses(info, info.free)
end
function freeparameterstuple(model::Type{<:CompositeModel})
    infos = getinfo(model)
    reduce(vcat, map(i -> _parameter_lenses(i, i.free), infos))
end

function frozenparameterstuple(model::Type{<:AbstractSpectralModel})
    info = getinfo(model)
    _parameter_lenses(info, info.frozen)
end
function frozenparameterstuple(model::Type{<:CompositeModel})
    infos = getinfo(model)
    reduce(vcat, map(i -> _parameter_lenses(i, i.frozen), infos))
end

function _parameter_lenses(info::ModelInfo, symbols)
    map(symbols) do s
        :(getproperty($(info.lens), $(Meta.quot(s))))
    end
end

end # module
