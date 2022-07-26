struct ProcessedSpectralModel{V}
    expression::Expr
    parameters::V
    multiplicatives::Vector{Expr}
    additives::Vector{Expr}
    convolutionals::Vector{Expr}
end

function Base.show(io::IO, m::ProcessedSpectralModel)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == ProcessedSpectralModel
        print(io, "ProcessedSpectralModel(M=$(length(m.multiplicatives)),A=$(length(m.additives)),C=$(length(m.convolutionals)))")
    else
        print(io, "ProcessedSpectralModel:\n   $(m.expression)\n")
        println(io, "  Models:")
        for m in m.multiplicatives
            print(io, "   - $m\n")
        end
        for m in m.additives
            print(io, "   - $m\n")
        end
        for m in m.convolutionals
            print(io, "   - $m\n")
        end
        println(io, "  Parameters:")
        for p in m.parameters
            print(io, "   - $p\n")
        end
    end
end

function assemblefunc!(tracker, m::AbstractSpectralModel{Additive})
    push!(tracker.additive, m)
    Symbol('a', length(tracker.additive))
end
function assemblefunc!(tracker, m::AbstractSpectralModel{Multiplicative})
    push!(tracker.multiplicative, m)
    Symbol('m', length(tracker.multiplicative))
end
function assemblefunc!(tracker, m::AbstractSpectralModel{Convolutional})
    push!(tracker.convolutional, m)
    Symbol('c', length(tracker.convolutional))
end

function assemblefunc!(tracker, cm::CompositeSpectralModel{K1,K2}) where {K1,K2}
    # Expr(
    #     :call,
    #     :(.*),
    #     assemblefunc!(tracker, cm.m1),
    #     assemblefunc!(tracker, cm.m2)
    # )
    a = assemblefunc!(tracker, cm.m1)
    b = assemblefunc!(tracker, cm.m2)
    :($a .* $b)
end

function assemblefunc!(tracker, cm::CompositeSpectralModel{Additive,Additive})
    a = assemblefunc!(tracker, cm.m1)
    b = assemblefunc!(tracker, cm.m2)
    :($a .+ $b)
end

function addparams!(all_params, m::AbstractSpectralModel)
    map(fieldnames(typeof(m))) do p
        i = 1
        symb = Symbol(p, i)
        while symb in first.(all_params)
            i += 1
            symb = Symbol(p, i)
        end
        push!(all_params, symb => getproperty(m, p))
        symb
    end
end

function parse_models_with_params!(all_params, models, root_symb::Char)
    map(enumerate(models)) do (i, m)
        params = addparams!(all_params, m)
        instance = Expr(:call, typeof(m), params...)
        Expr(:call, :(=), Symbol(root_symb, i), instance)
    end
end

function processmodel(cm::AbstractSpectralModel)
    # assemble evaluation
    tracker = (
        additive = AbstractSpectralModel{Additive}[],
        multiplicative = AbstractSpectralModel{Multiplicative}[],
        convolutional = AbstractSpectralModel{Convolutional}[]
    )

    expr = assemblefunc!(tracker, cm)
    all_params = Pair{Symbol, Float64}[]

    mults = parse_models_with_params!(all_params, tracker.multiplicative, 'm')
    adds = parse_models_with_params!(all_params, tracker.additive, 'a')
    convs = parse_models_with_params!(all_params, tracker.convolutional, 'c')

    ProcessedSpectralModel(
        expr, all_params, mults, adds, convs
    )
end

export processmodel, ProcessedSpectralModel
