UNTRACKED_ERROR = zeros(Float64, 1)

function resize_untracked_error(n)
    resize!(UNTRACKED_ERROR, n)
end

macro wrap_xspec_model_ccall(
    func_name,
    input,
    output,
    err,
    parameters,
    spec_number,
    init_string,
)
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
    end |> esc
end

macro xspecmodel(model_kind, func_name, model)
    model_args = model.args[3]
    model_name = model.args[2].args[1]
    model_type_params = model.args[2].args[2]
    model_args_symbols = [
        :($(i.args[1].args[1])) for i in model_args.args if (i isa Expr) && (i.head == :(=))
    ]
    # remove normalisation as a parameter from Additive models
    if model_kind == :Additive
        deleteat!(model_args_symbols, 1)
    end

    parsed_model_args = [:(get_value($i)) for i in model_args_symbols]

    quote
        @with_kw struct $(model_name){$(model_type_params)} <: AbstractSpectralModel
            $(model.args[3].args...)
        end

        modelkind(::Type{<:$(model_name)}) = $(model_kind)()
        implementation(::Type{<:$(model_name)}) = XSPECImplementation()

        function invoke!(
            flux::AbstractArray,
            energy::AbstractArray,
            m::Type{<:$(model_name)},
            $(model_args_symbols...);
            spectral_number = 1,
            init_string = "",
        )
            @assert length(flux) + 1 == length(energy)
            if length(UNTRACKED_ERROR) < length(flux)
                resize_untracked_error(length(flux))
            end

            @wrap_xspec_model_ccall(
                $(func_name),
                energy,
                flux,
                UNTRACKED_ERROR,
                [$(parsed_model_args...)],
                spectral_number,
                init_string
            )
        end
    end |> esc
end
