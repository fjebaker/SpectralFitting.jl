
macro wrap_xspec_model_ccall(func_name, input, output, err, parameters, spec_number, init_string)
    quote
        ccall(
            ($(func_name), libXSFunctions),
            Cvoid,
            (Ref{Float64}, Int32, Ref{Float64}, Int32, Ref{Float64}, Ref{Float64}, Cstring),
            $(input),
            length($(input)) - 1,
            $(parameters),
            $(spec_number),
            $(output),
            $(err),
            $(init_string),
        )
        # double up the last value
        # since XSPEC expecty length(energy) == length(flux) + 1
        $(output)[end] = $(output)[end-1]
        $(output)
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

    # remove normalisation as a parameter from Additive models
    if model_kind == :Additive
        deleteat!(model_args_symbols, 1)
    end

    quote
        @with_kw struct $(model_name) <: AbstractSpectralModel{$(model_kind)}
            $(model.args[3].args...)
        end

        function invokemodel(
            input::AbstractArray,
            m::$(model_name);
            spectral_number = 1,
            init_string = "",
        )
            output = zeros(eltype(input), size(input))
            err = zeros(eltype(input), size(input))

            @wrap_xspec_model_ccall(
                $(func_name),
                input,
                output,
                err,
                [$(model_args_symbols...)],
                spectral_number,
                init_string
            )
        end

        function invokemodel!(
            output::AbstractArray,
            err::AbstractArray,
            input::AbstractArray,
            m::$(model_name);
            spectral_number = 1,
            init_string = "",
        )
            @assert size(output) == size(input) "Output and input must have same dimensions."
            @assert size(err) == size(input) "Error array and input must have same dimensions."

            @wrap_xspec_model_ccall(
                $(func_name),
                input,
                output,
                err,
                [$(model_args_symbols...)],
                spectral_number,
                init_string
            )
        end
    end |> esc
end
