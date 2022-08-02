
function add_flux_count!(flux_count, symb)
    flux_count[] += 1
    flux = Symbol(:flux, flux_count[])
    :(invokemodel!($flux, energy, $symb))
end

function resolve_flux_combine!(flux_count, op)
    flux_right = Symbol(:flux, flux_count[])
    flux_left = Symbol(:flux, flux_count[] - 1)
    flux_count[] -= 1
    expr = Expr(:call, op, flux_left, flux_right)
    :(@.($flux_left = $expr))
end

function __build_statements(expression, models)
    flux_count = Ref(0)
    statements = Expr[]
    recursive_expression_call(expression) do (left, right, op)
        @assert flux_count[] < 4 error("Flux count exceeded 3.")

        if typeof(right) === Symbol
            M = modelkind(models[right])
            @assert (M === Additive) error("Right, if symbol, should be Additive.")
            push!(statements, add_flux_count!(flux_count, right))
        end

        if typeof(left) === Symbol
            M = modelkind(models[left])
            if M === Multiplicative
                push!(statements, add_flux_count!(flux_count, left))
                push!(statements, resolve_flux_combine!(flux_count, op))
            elseif M === Convolutional
                flux = Symbol(:flux, flux_count[])
                push!(statements, :(invokemodel!($flux, energy, $left)))
            else
                error("Unrecognised model type $M.")
            end

        elseif isnothing(left)
            push!(statements, resolve_flux_combine!(flux_count, :(+)))
        else
            error("Left should always be a symbol or nothing.")
        end

        nothing
    end

    statements
end

function __build_statements(symbol::Symbol, models)
    [:(invokemodel!(flux1, energy, $symbol))]
end

function __build_model_instance(models, modelparams)
    map(models) do (s, M)
        params = modelparams[s]
        model_name = Base.typename(M).name
        :($s = $(model_name)($(params...)))
    end
end

function __build_parameter_statements(params)
    fps = filter(!isfrozen, last.(params))

    i = 0
    statements = map(params) do (s, p)
        if p.frozen
            val = value(p)
            :($s = $val)
        else
            i += 1
            :($s = params[$i])
        end
    end

    statements, fps
end

function build_simple_no_eval(psm::ProcessedSpectralModel)
    statements = __build_statements(psm.expression, Dict(psm.models))
    model_instances = __build_model_instance(psm.models, psm.modelparams)
    parameters, fps = __build_parameter_statements(psm.parameters)

    func = :((flux1, flux2, flux3, energy, params) -> begin
        @fastmath @inbounds begin
            $(parameters...)
            $(model_instances...)
            $(statements...)
            return flux1
        end
    end)

    func, fps
end

function build_simple(psm::ProcessedSpectralModel)
    func, fps = build_simple_no_eval(psm)
    evaled_func = eval(func)
    closure_wrapper = (args...) -> Base.invokelatest(evaled_func, args...)
    closure_wrapper, fps
end

function __build_parameter_distribution_statements(params)
    i = 0
    statements = map(params) do (s, p)
        if p.frozen
            val = value(p)
            :($s = $val)
        else
            i += 1
            distr = if typeof(p) <: FitParameter
                # use default distribution
                val = value(p)
                lb = lowerbound(p)
                ub = upperbound(p) == Inf ? 1e6 : upperbound(p)
                :(Turing.TruncatedNormal($val, 2.0, $lb, $ub))
            elseif typeof(p) <: AbstractFitDistributionParameter
                # ... TODO
                p
            else
                error("Unknown fit parameter type $(typeof(p)).")
            end
            :($s ~ $distr)
        end
    end
    statements
end

function build_turing_no_eval(psm::ProcessedSpectralModel)
    statements = __build_statements(psm.expression, Dict(psm.models))
    model_instances = __build_model_instance(psm.models, psm.modelparams)
    parameters = __build_parameter_distribution_statements(psm.parameters)

    :(
        begin
            (target, energy, error, channels, flux1, flux2, flux3) -> begin
                Turing.@model function __wrapped_mcmc_func(target)
                    $(parameters...)
                    $(model_instances...)
                    $(statements...)

                    for (i, j) in enumerate(channels)
                        target[i] ~ Turing.Normal(flux1[j], error[i])
                    end
                end

                __wrapped_mcmc_func(target)
            end
        end
    )
end

function build_turing(psm::ProcessedSpectralModel)
    evaled_func = build_turing_no_eval(psm) |> eval
    closure_wrapper = (args...) -> Base.invokelatest(evaled_func, args...)
end

export build_simple, build_simple_no_eval, build_turing, build_turing_no_eval
