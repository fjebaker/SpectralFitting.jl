
macro wrap_xspec_model_ccall(func_name, input, parameters, spec_number, init_string)
    quote
        output = zeros(eltype($(input)), size($(input)))
        err = zeros(eltype($(input)), size($(input)))
        ccall(
            ($(func_name), libXSFunctions),
            Cvoid,
            (Ref{Float64}, Int32, Ref{Float64}, Int32, Ref{Float64}, Ref{Float64}, Cstring),
            $(input),
            length($(input)) - 1,
            $(parameters),
            $(spec_number),
            output,
            err,
            $(init_string),
        )
        # double up the last value
        # since XSPEC expecty length(energy) == length(flux) + 1
        output[end] = output[end-1]
        output
    end |> esc
end

macro xspecmodel(func_name, model)
    model_args = model.args[3]
    model_name = model.args[2].args[1]
    model_kind = model.args[2].args[2]
    model_args_symbols = [
        :(m.$(i.args[1].args[1])) for
        i in model_args.args if (i isa Expr) && (i.head == :(=))
    ]

    if model_kind == :Additive
        deleteat!(model_args_symbols, 1)
    end

    quote
        @with_kw struct $(model_name) <: AbstractSpectralModel{$(model_kind)}
            $(model.args[3].args...)
        end

        function invokemodel(
            m::$(model_name),
            input::AbstractArray;
            spectral_number = 1,
            init_string = "",
        )
            @wrap_xspec_model_ccall(
                $(func_name),
                input,
                [$(model_args_symbols...)],
                spectral_number,
                init_string
            )
        end
    end |> esc
end
