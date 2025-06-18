export FitResult, update_model!, get_objective, get_objective_variance

struct FitResult{Config<:FittingConfig,U,Err,T,Sol}
    config::Config
    u::U
    err::Err
    stats::Vector{T}
    sol::Sol
end

result_count(r::FitResult) = length(r.config.prob.model.m)

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
    model = get_model(slice)
    print(io, "Model: ")
    printstyled(io, _model_name(model), color = :cyan)
    println(io)
    print(io, " . Name : ")

    params, syms = _all_parameters_with_names(model)
    free_syms = syms[isfree.(params)]

    param_padding = max(10, maximum(length, free_syms)) + 1

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
    bindings = result.config.parameter_bindings[i]
    mask = result.config.parameter_cache.free_mask[bindings]

    err_slice = if isnothing(result.err)
        nothing
    else
        err_parameters = update_free_parameters!(result.config.parameter_cache, result.err)
        err_parameters[bindings][mask]
    end

    all_parameters = update_free_parameters!(result.config.parameter_cache, result.u)
    u_slice = all_parameters[bindings][mask]

    FitResultSlice(i, result, u_slice, err_slice, result.stats[i])
end

function calculate_objective!(slice::FitResultSlice, u0)
    I = slice.parent.config.parameter_bindings[slice.index]
    mask = slice.parent.config.parameter_cache.free_mask[I]
    @assert count(mask) == length(u0)

    all_parameters = _get_parameters(slice.parent.config.parameter_cache, u0)
    # update the free parameters
    @views all_parameters[I][mask] .= u0

    calculate_objective!(slice.parent.config, all_parameters, slice.index)
end

_get_data_cache(slice::FitResultSlice) = slice.parent.config.data_cache[slice.index]
get_model(slice::FitResultSlice) = slice.parent.config.prob.model.m[slice.index]
get_dataset(slice::FitResultSlice) = slice.parent.config.prob.data.d[slice.index]
get_objective(slice::FitResultSlice) = _get_data_cache(slice).objective
get_objective_variance(slice::FitResultSlice) = _get_data_cache(slice).variance
plotting_domain(slice::FitResultSlice) = plotting_domain(get_dataset(slice))

function finalize_result(config::FittingConfig, params, sol; σparams = nothing)
    I = ((1:model_count(config.prob))...,)

    measures = map(I) do i
        obj = measure_objective!(config, params, i)
    end |> collect

    FitResult(config, params, σparams, measures, sol)
end

function measure(s::AbstractStatistic, slice::FitResultSlice, u = slice.u)
    ŷ = calculate_objective!(slice, u)
    measure(s, get_objective(slice), ŷ, get_objective_variance(slice))
end

function measure(stat::AbstractStatistic, result::FitResult, args...; kwargs...)
    measure(stat, result[1], args...; kwargs...)
end

update_model!(model::AbstractSpectralModel, result::FitResult) =
    update_model!(model, result[1])

function residuals(slice::FitResultSlice)
    y = calculate_objective!(slice, slice.u)
    obj, var = get_objective(slice), get_objective_variance(slice)
    @. (obj - y) / sqrt(var)
end

function update_model!(model::AbstractSpectralModel, result::FitResultSlice)
    all_params = update_free_parameters!(result.parent.config.parameter_cache, result.u)
    for (p, r) in zip(parameter_vector(model), result.u)
        set_value!(p, r)
    end
    model
end
