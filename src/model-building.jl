
function resolve_flux_combine!(flux_count, op)
    i = (flux_count[] -= 1)
    flux_right = Symbol(:flux, i + 1)
    flux_left = Symbol(:flux, i)
    expr = Expr(:call, op, flux_left, flux_right)
    :(@.($flux_left = $expr))
end

function assemble_invoke(flux, symb, M, params)
    base_name = model_type_name(M)
    :(invokemodel!($flux, energy, $base_name{Float64}, $(params...)))
end

function add_flux_invokation(flux_count, symb, params, M, ::Convolutional)
    # does not increase flux count
    flux = Symbol(:flux, flux_count[])
    assemble_invoke(flux, symb, M, params)
end

function add_flux_invokation(flux_count, symb, params, M, ::AbstractSpectralModelKind)
    i = (flux_count[] += 1)
    flux = Symbol(:flux, i)
    assemble_invoke(flux, symb, M, params)
end

function assemble_invokation_statements(expression, s_to_M, s_to_params)
    flux_count = Ref(0)
    statements = Expr[]
    recursive_expression_call(expression) do (left, right, op)
        @assert (0 ≤ (flux_count[]) < 4) error("Flux counts out of bounds $(flux_count[]).")
        if right isa Symbol
            M = s_to_M[right]
            K = modelkind(M)
            push!(
                statements,
                add_flux_invokation(flux_count, right, s_to_params[right], M, K),
            )
        end
        if left isa Symbol
            M = s_to_M[left]
            K = modelkind(M)
            push!(
                statements,
                add_flux_invokation(flux_count, left, s_to_params[left], M, K),
            )
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
function assemble_invokation_statements(symb::Symbol, s_to_M, s_to_params)
    [assemble_invoke(:flux1, symb, s_to_M[symb], s_to_params[symb])]
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

function assemble_instance_arguments(instances)
    additional_args = Symbol[]
    for (s, m) in instances
        push!(additional_args, make_additional_invoke_parameters_symbols(m, s)...)
    end
    additional_args
end

function build_simple_no_eval(psm::ProcessedSpectralModel)
    statements = assemble_invokation_statements(
        psm.expression,
        Dict(psm.model_type_index),
        psm.modelparams,
    )
    parameters, fps = assemble_parameter_assignments(psm.parameters)
    additional_args = assemble_instance_arguments(psm.instances)
    func = :(
        (flux1, flux2, flux3, energy, params, $(additional_args...)) -> begin
            @fastmath @inbounds begin
                $(parameters...)
                $(statements...)
                return flux1
            end
        end
    )

    Base.remove_linenums!(func)

    func, fps
end

function assemble_additional_args(instances)
    args = []
    for (_, i) in instances
        push!(args, additional_invoke_parameters(i)...)
    end
    args
end

function runtime_eval_and_wrap(expr, instances)
    f = RuntimeGeneratedFunctions.@RuntimeGeneratedFunction(expr)
    additional_args = assemble_additional_args(instances)
    if !isempty(additional_args)
        (args...) -> begin
            f(args..., additional_args...)
        end
    else
        f
    end
end

function build_simple(psm::ProcessedSpectralModel)
    func, fps = build_simple_no_eval(psm)
    runtime_eval_and_wrap(func, psm.instances), fps
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
    statements = assemble_invokation_statements(
        psm.expression,
        Dict(psm.model_type_index),
        psm.modelparams,
    )
    parameters = assemble_parameters_as_distributions(psm.parameters)
    additional_args = assemble_instance_arguments(psm.instances)
    func = :(
        begin
            (energy, target, error, channels, flux1, flux2, flux3, $(additional_args...)) -> begin
                Turing.@model function __wrapped_mcmc_func(target)
                    $(parameters...)
                    $(statements...)

                    for (i, j) in enumerate(channels)
                        target[i] ~ Turing.Normal(flux1[j], error[i])
                    end
                end

                __wrapped_mcmc_func(target)
            end
        end
    )

    Base.remove_linenums!(func)

    func
end

function build_turing(psm::ProcessedSpectralModel)
    runtime_eval_and_wrap(build_turing_no_eval(psm), psm.instances)
end

export build_simple, build_simple_no_eval, build_turing, build_turing_no_eval
