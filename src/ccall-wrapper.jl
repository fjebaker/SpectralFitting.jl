UNTRACKED_ERROR = zeros(Float64, 1)

function resize_untracked_error(n)
    resize!(UNTRACKED_ERROR, n)
end

_unsafe_unwrap_parameter(ptr, N, ::AbstractSpectralModelKind, T::Type) = 
    unsafe_wrap(Vector{T}, ptr, N, own = false)
# additive models don't include normalisation, so we offset
_unsafe_unwrap_parameter(ptr, N, ::Additive, T::Type) = 
    unsafe_wrap(Vector{T}, ptr + sizeof(T), N - 1, own = false)

@inline function unsafe_parameter_vector(
    alloc_model::Vector{<:AbstractSpectralModel{T}},
) where {T}
    model = first(alloc_model)
    N = length(fieldnames(typeof(model)))

    model_ptr = pointer(alloc_model)
    ptr = Base.unsafe_convert(Ptr{T}, model_ptr)
    # reinterpret as an unowned vector of `T`
    _unsafe_unwrap_parameter(ptr, N, modelkind(model), T)
end

# pointer hackery
function unsafe_parameter_vector_conditioned(
    alloc_model::Vector{<:AbstractSpectralModel{T}},
) where {T}
    model = first(alloc_model)
    if implementation(typeof(model)) == JuliaImplementation()
        throw("This method is only for `XSPECImplementation` models.")
    end
    if !isbitstype(T)
        throw("Can only reinterpret to bits type (convert `$T`).")
    end
    unsafe_parameter_vector(alloc_model)
end

macro wrap_xspec_model_ccall(
    func_name,
    callsite,
    input,
    output,
    err,
    parameters,
    spec_number,
    init_string,
)
    quote
        ccall(
            ($(func_name), $(callsite)),
            Cvoid,
            (Ptr{Float64}, Cint, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Float64}, Cstring),
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
    @xspecmodel model_kind func_name model

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

The macro will then generate the following functions
```julia
modelkind(::Type{<:XS_PowerLaw})
implementation(::Type{<:XS_PowerLaw})
invoke!(::Type{<:XS_PowerLaw})
```
"""
macro xspecmodel(c_function, model)
    model_args = model.args[3]
    model_name = model.args[2].args[1].args[1]
    model_type_params = model.args[2].args[1].args[2:end]
    model_kind = model.args[2].args[end]
    symbols = [i.args[1] for i in model_args.args if (i isa Expr) && (i.head == :(::))]

    if model_kind.args[2] == :Additive
        # get rid of the normalisation from function arguments
        symbols = symbols[2:end]
    end

    if c_function isa QuoteNode
        func_name = c_function
        callsite = libXSFunctions
    else
        func_name = c_function.args[1]
        callsite = c_function.args[2]
    end


    quote
        Base.@__doc__ struct $(model_name){$(model_type_params...)} <: $(model_kind)
            $(model_args.args...)
        end

        SpectralFitting.implementation(::Type{<:$(model_name)}) = XSPECImplementation()

        function SpectralFitting.invoke!(
            flux::AbstractArray,
            energy::AbstractArray,
            model::M,
            spectral_number = 1,
            init_string = "",
        ) where {M<:$(model_name)}
            @assert length(flux) + 1 == length(energy)
            SpectralFitting.ensure_model_data(M)

            if length(SpectralFitting.UNTRACKED_ERROR) < length(flux)
                SpectralFitting.resize_untracked_error(length(flux))
            end

            # allocate the model
            alloc_model = [model]
            params = SpectralFitting.unsafe_parameter_vector_conditioned(alloc_model)

            @wrap_xspec_model_ccall(
                $(func_name),
                $(callsite),
                energy,
                flux,
                SpectralFitting.UNTRACKED_ERROR,
                params,
                spectral_number,
                init_string
            )
        end
        $(model_name)
    end |> esc
end

export @xspecmodel, @wrap_xspec_model_ccall
