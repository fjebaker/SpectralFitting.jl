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

# include utility functions after struct definition
include("parsing-utilities.jl")

function new_closure_param!(ga::GenerationAggregate, p::Symbol)
    param = Base.gensym(p)
    push_closure_param!(ga, param)
    param
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
            @inbounds let ($(flux_unpack...))
                $(closures...)
                $(p_assignments...)
                $(statements...)
                return flux1
            end
        end
    end
end

function make_flux_unpack(N)
    symbols = (Symbol(:flux, i) for i = 1:N)
    unpacks = [:($f = view(fluxes, :, $i)) for (i, f) in enumerate(symbols)]
    unpacks
end

function generated_model_call!(fluxes, energy, model, free_params, frozen_params)
    # propagate information about free parameters to allow for AD
    ga = assemble_aggregate_info(model, eltype(free_params))
    flux_unpack = make_flux_unpack(ga.maximum_flux_count)
    p_assign = assemble_parameter_assignment(ga, model)
    closures = assemble_closures(ga, model)
    generate_call(flux_unpack, closures, p_assign, ga.statements)
end
function generated_model_call!(fluxes, energy, model, params)
    # propagate information about free parameters to allow for AD
    ga = assemble_aggregate_info(model, eltype(params))
    flux_unpack = make_flux_unpack(ga.maximum_flux_count)
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


function assemble_aggregate_info(model::Type{<:AbstractSpectralModel}, NumType)
    ga = FunctionGeneration.GenerationAggregate(NumType)
    _addinfoinvoke!(ga, model, :(model))
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

function model_parameters_tuple(model::Type{<:AbstractSpectralModel})
    info = getinfo(model)
    _parameter_lens(info, info.symbols)
end
function model_parameters_tuple(model::Type{<:CompositeModel})
    infos = getinfo(model)
    reduce(vcat, map(i -> _parameter_lens(i, i.symbols), infos))
end

function free_parameters_tuple(model::Type{<:AbstractSpectralModel})
    info = getinfo(model)
    _parameter_lens(info, info.free)
end
function free_parameters_tuple(model::Type{<:CompositeModel})
    infos = getinfo(model)
    reduce(vcat, map(i -> _parameter_lens(i, i.free), infos))
end

function frozen_parameters_tuple(model::Type{<:AbstractSpectralModel})
    info = getinfo(model)
    _parameter_lens(info, info.frozen)
end
function frozen_parameters_tuple(model::Type{<:CompositeModel})
    infos = getinfo(model)
    reduce(vcat, map(i -> _parameter_lens(i, i.frozen), infos))
end

function _parameter_lens(info::ModelInfo, symbols)
    map(symbols) do s
        :(getproperty($(info.lens), $(Meta.quot(s))))
    end
end

end # module
