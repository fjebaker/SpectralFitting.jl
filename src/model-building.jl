function add_calculate_flux!(flux_count, symb, ::Convolutional)
    # does not increase flux count
    flux = Symbol(:flux, flux_count[])
    :(invokemodel!($flux, energy, $symb))
end

function add_calculate_flux!(flux_count, symb, ::AbstractSpectralModelKind)
    i = (flux_count[] += 1)
    flux = Symbol(:flux, i)
    :(invokemodel!($flux, energy, $symb))
end

function resolve_flux_combine!(flux_count, op)
    i = (flux_count[] -= 1)
    flux_right = Symbol(:flux, i + 1)
    flux_left = Symbol(:flux, i)
    expr = Expr(:call, op, flux_left, flux_right)
    :(@.($flux_left = $expr))
end

function assemble_execution_statements(expression, models)
    flux_count = Ref(0)
    statements = Expr[]
    recursive_expression_call(expression) do (left, right, op)
        @assert (0 ≤ (flux_count[]) < 4) error("Flux counts out of bounds $(flux_count[]).")
        if right isa Symbol
            M = modelkind(models[right])
            push!(statements, add_calculate_flux!(flux_count, right, M))
        end
        if left isa Symbol
            M = modelkind(models[left])
            push!(statements, add_calculate_flux!(flux_count, left, M))
        end
        # catch case for Convolution operations
        if !isnothing(op)
            push!(statements, resolve_flux_combine!(flux_count, op))
        end
        nothing
    end
    statements
end


# edge-case when only one model
function assemble_execution_statements(symbol::Symbol, _)
    [:(invokemodel!(flux1, energy, $symbol))]
end

function assemble_model_constructors(models, modelparams)
    map(models) do (s, M)
        params = modelparams[s]
        model_name = model_type_name(M)
        :($s = $(model_name)($(params...)))
    end
end

function assemble_parameter_assignments(params)
    i = 0
    statements = map(params) do (s, p)
        if is_frozen(p)
            val = get_value(p)
            :($s = $val)
        else
            i += 1
            :($s = params[$i])
        end
    end
    statements, filter(!is_frozen, last.(params))
end

function build_simple_no_eval(psm::ProcessedSpectralModel)
    statements = assemble_execution_statements(psm.expression, Dict(psm.models))
    model_instances = assemble_model_constructors(psm.models, psm.modelparams)
    parameters, fps = assemble_parameter_assignments(psm.parameters)

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

function eval_and_wrap(expr)
    f = eval(expr)
    (args...; kwargs...) -> Base.invokelatest(f, args...; kwargs...)
end

function build_simple(psm::ProcessedSpectralModel)
    func, fps = build_simple_no_eval(psm)
    eval_and_wrap(func), fps
end

function assemble_param_as_distribution(p::FitParam)
    type, args = as_distribution(p)
    :($(type)($(args...)))
end

function assemble_parameters_as_distributions(params)
    i = 0
    statements = map(params) do (s, p)
        if is_frozen(p)
            val = get_value(p)
            :($s = $val)
        else
            i += 1
            distr = assemble_param_as_distribution(p)
            :($s ~ $distr)
        end
    end
    statements
end

function build_turing_no_eval(psm::ProcessedSpectralModel)
    statements = assemble_execution_statements(psm.expression, Dict(psm.models))
    model_instances = assemble_model_constructors(psm.models, psm.modelparams)
    parameters = assemble_parameters_as_distributions(psm.parameters)
    :(
        begin
            (energy, target, error, channels, flux1, flux2, flux3) -> begin
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
    eval_and_wrap(build_turing_no_eval(psm))
end

export build_simple, build_simple_no_eval, build_turing, build_turing_no_eval
