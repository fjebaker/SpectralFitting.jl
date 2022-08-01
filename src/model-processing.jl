mutable struct FitParameter{T}
    val::T
    lower_bound::T
    upper_bound::T
    frozen::Bool

    FitParameter(val::T; lower_bound::T = T(0.0), upper_bound::T = T(Inf), frozen::Bool = false) where {T} = new{T}(val, lower_bound, upper_bound, frozen)
end

value(f::FitParameter) = f.val
lowerbound(f::FitParameter) = f.lower_bound
upperbound(f::FitParameter) = f.upper_bound

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, m::FitParameter{T}) where {T}
    Base.show(io, m)
end

function Base.show(io::IO, m::FitParameter{T}) where {T}
    s = "$(m.val) (lb: $(m.lower_bound), ub: $(m.upper_bound))"
    if m.frozen
        s *= " [frozen]"
    end
    print(io, s)
end

struct ProcessedSpectralModel{V}
    expression::Expr
    parameters::V
    multiplicatives::Vector{Pair{Symbol,Expr}}
    additives::Vector{Pair{Symbol,Expr}}
    convolutionals::Vector{Pair{Symbol,Expr}}
end

function getmodels(pm::ProcessedSpectralModel)
    Iterators.flatten((pm.multiplicatives, pm.additives, pm.convolutionals))
end

function getparams(pm::ProcessedSpectralModel)
    pm.parameters
end

function getfreeparams(pm::ProcessedSpectralModel)
    filter(p -> !last(p).frozen, pm.parameters)
end

function getfrozenparams(pm::ProcessedSpectralModel)
    filter(p -> last(p).frozen, pm.parameters)
end

function freeze!(pm::ProcessedSpectralModel, p::Symbol)
    index = findfirst(i -> first(i) == p, pm.parameters)
    last(pm.parameters[index]).frozen = true
end

freeze!(pm::ProcessedSpectralModel, ps::Symbol...) = foreach(p -> freeze!(pm, p), ps)

function unfreeze!(pm::ProcessedSpectralModel, p::Symbol)
    index = findfirst(i -> first(i) == p, pm.parameters)
    last(pm.parameters[index]).frozen = false
end

unfreeze!(pm::ProcessedSpectralModel, ps::Symbol...) = foreach(p -> unfreeze!(pm, p), ps)


function Base.show(io::IO, m::ProcessedSpectralModel)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == ProcessedSpectralModel
        print(
            io,
            "ProcessedSpectralModel(M=$(length(m.multiplicatives)),A=$(length(m.additives)),C=$(length(m.convolutionals)))",
        )
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
    :($a * $b)
end

function assemblefunc!(tracker, cm::CompositeSpectralModel{Additive,Additive})
    a = assemblefunc!(tracker, cm.m1)
    b = assemblefunc!(tracker, cm.m2)
    :($a + $b)
end

function addparams!(all_params, m::AbstractSpectralModel)
    map(fieldnames(typeof(m))) do p
        i = 1
        symb = Symbol(p, '_', i)
        while symb in first.(all_params)
            i += 1
            symb = Symbol(p, '_', i)
        end
        push!(all_params, symb => FitParameter(getproperty(m, p)))
        symb
    end
end

function parse_models_with_params!(all_params, models, root_symb::Char)
    map(enumerate(models)) do (i, m)
        params = addparams!(all_params, m)
        instance = Expr(:call, typeof(m), params...)
        Symbol(root_symb, i) => instance
    end
end

function processmodel(cm::AbstractSpectralModel)
    # assemble evaluation
    tracker = (
        additive = AbstractSpectralModel{Additive}[],
        multiplicative = AbstractSpectralModel{Multiplicative}[],
        convolutional = AbstractSpectralModel{Convolutional}[],
    )

    expr = assemblefunc!(tracker, cm)
    all_params = Pair{Symbol,FitParameter{Float64}}[]

    mults = parse_models_with_params!(all_params, tracker.multiplicative, 'm')
    adds = parse_models_with_params!(all_params, tracker.additive, 'a')
    convs = parse_models_with_params!(all_params, tracker.convolutional, 'c')

    ProcessedSpectralModel(expr, all_params, mults, adds, convs)
end

function build_lsqfit(pm::ProcessedSpectralModel)
    flux_calls = [:($name = $(model)(energy)) for (name, model) in getmodels(pm)]

    i = 1
    params = map(pm.parameters) do (p, f)
        if f.frozen
            val = value(f)
            :($p = $val)
        else
            i += 1
            :($p = params[$(i-1)])
        end
    end

    # params = [:($p = params[$i]) for (i, p) in enumerate(first.(pm.parameters))]
    func = :((energy, params) -> begin
        $(params...)
        $(flux_calls...)
        res = @. $(pm.expression)
        @view res[1:end-1]
    end) |> eval

    fs = last.(getfreeparams(pm))
    # model, p0, lower bounds, upper bounds
    func, value.(fs), lowerbound.(fs), upperbound.(fs)
end

export processmodel, ProcessedSpectralModel, getmodels, build_lsqfit, FitParameter, value, freeze!, unfreeze!
