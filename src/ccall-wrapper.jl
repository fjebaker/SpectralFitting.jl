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

"""
    @xspecmodel model_kind, func_name, model

Used to wrap additional XSPEC models, generating the needed [`AbstractSpectralModel`](@ref)
implementation.

# Examples

```julia
@xspecmodel Additive :C_powerlaw struct XS_PowerLaw{F1,F2}
    "Normalisation."
    K::F1 = FitParam(1.0)
    "Photon index."
    a::F2 = FitParam(0.5)
end
```

We define a new structure `XS_PowerLaw` with two parameters, but since the model is [`Additive`](@ref),
only a single parameter (`a`) is passed to the XSPEC function. The function we bind to this model
is `:C_powerlaw` from the XSPEC C wrappers.

The macro will, in this case, generate the following functions
```julia
modelkind(::Type{<:XS_PowerLaw})
implementation(::Type{<:XS_PowerLaw})
invoke!(::Type{<:XS_PowerLaw})
```
"""
macro xspecmodel(model_kind, func_name, model)
    model_args = model.args[3]
    model_name = model.args[2].args[1]
    model_type_params = model.args[2].args[2:end]
    model_args_symbols = [
        :($(i.args[1].args[1])) for i in model_args.args if (i isa Expr) && (i.head == :(=))
    ]
    # remove normalisation as a parameter from Additive models
    if model_kind == :Additive
        deleteat!(model_args_symbols, 1)
    end

    parsed_model_args = [:(get_value($i)) for i in model_args_symbols]

    quote
        @with_kw struct $(model_name){$(model_type_params...)} <: AbstractSpectralModel
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


export @xspecmodel
