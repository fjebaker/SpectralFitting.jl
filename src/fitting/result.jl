export FitResult

struct FitResult{Config<:FittingConfig,U,Err,T,Sol}
    config::Config
    u::U
    err::Err
    stats::Vector{T}
    sol::Sol
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(result::FitResult))
    buff = IOBuffer()
    buff_c = IOContext(buff, io)

    total_stat = prettyfloat(sum(result.stats))

    println(buff_c, "FitResult:")
    print(buff_c, " ")
    for i = 1:length(result.stats)
        slice = result[i]

        buff2 = IOBuffer()
        buff2_c = IOContext(buff2, io)
        _pretty_print_result(buff2_c, slice)

        b = String(take!(buff2))
        r = indent(b, 1)

        print(buff_c, r)
    end
    text = String(take!(buff))
    print(
        io,
        encapsulate(text) *
        "Σ$(statistic_symbol(fit_statistic(result.config))) = $(total_stat)",
    )
end

calculate_objective!(result::FitResult, u0) = calculate_objective!(result.config, u0)

struct FitResultSlice{P<:FitResult,U,Err,T}
    index::Int
    parent::P
    u::U
    err::Err
    stats::T
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(slice::FitResultSlice))
    buff = IOBuffer()
    buff_c = IOContext(buff, io)

    println(buff_c, "FitResultSlice:")
    _pretty_print_result(buff_c, slice)

    text = String(take!(buff))
    print(io, encapsulate(text))
end

function _pretty_print_result(io::IO, slice::FitResultSlice)
    param_padding = 10
    model = get_model(slice)
    print(io, "Model: ")
    printstyled(io, _model_name(model), color = :cyan)
    println(io)
    print(io, " . Name : ")

    params, syms = _all_parameters_with_symbols(model)
    free_syms = syms[isfree.(params)]

    for s in free_syms
        print(io, rpad(s, param_padding))
    end
    println(io)

    print(io, " . u    : ")
    for v in slice.u
        print(io, rpad(prettyfloat(v), param_padding))
    end
    println(io)

    print(io, " . Δu   : ")
    if !isnothing(slice.err)
        for v in slice.err
            print(io, rpad(prettyfloat(v), param_padding))
        end
    else
        printstyled(io, "nothing", color = :gray)
    end
    println(io)

    stat_sym = statistic_symbol(fit_statistic(slice.parent.config))
    print(io, " . $(rpad(stat_sym, 6 - length(stat_sym))) : $(prettyfloat(slice.stats))")
    println(io)
end

function Base.getindex(result::FitResult, i)
    all_parameters = update_free_parameters!(result.config.parameter_cache, result.u)
    mask = result.config.parameter_cache.free_mask[result.config.parameter_bindings[i]]
    FitResultSlice(
        i,
        result,
        all_parameters[result.config.parameter_bindings[i]][mask],
        isnothing(result.err) ? nothing : result.err[result.config.bindings[i]],
        result.stats[i],
    )
end

function calculate_objective!(slice::FitResultSlice, u0)
    I = slice.parent.config.parameter_bindings[slice.index]
    mask = slice.parent.config.parameter_cache.free_mask[I]
    @assert count(mask) == length(u0)

    all_parameters = _get_parameters(slice.parent.config.parameter_cache, u0)
    # update the free parameters
    @views all_parameters[mask] .= u0

    calculate_objective!(slice.parent.config, all_parameters, slice.index)
end

_get_data_cache(slice::FitResultSlice) = slice.parent.config.data_cache[slice.index]
get_model(slice::FitResultSlice) = slice.parent.config.prob.model.m[slice.index]
get_objective(slice::FitResultSlice) = _get_data_cache(slice).objective
get_objective_variance(slice::FitResultSlice) = _get_data_cache(slice).variance

function finalize_result(config::FittingConfig, params, sol; σparams = nothing)
    I = ((1:model_count(config.prob))...,)

    measures = map(I) do i
        obj = measure_objective!(config, params, i)
    end |> collect

    FitResult(config, params, σparams, measures, sol)
end
