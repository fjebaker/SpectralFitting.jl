struct BinnedData{V} <: AbstractDataset
    domain::V
    codomain::V
    domain_variance::Union{Nothing,V}
    codomain_variance::Union{Nothing,V}
    name::String
end

function BinnedData(
    domain,
    codomain;
    domain_variance = nothing,
    codomain_variance = nothing,
    name = "[no-name]",
)
    @assert length(domain) == length(codomain) + 1 "Data does not look like it is binned!"
    BinnedData(domain, codomain, domain_variance, codomain_variance, name)
end

supports(::Type{<:BinnedData}) = (ContiguouslyBinned(), OneToOne())

bin_widths(dataset::BinnedData) = diff(dataset.domain)
objective_units(::BinnedData) = u"counts / (s * keV)"
spectrum_energy(dataset::BinnedData) =
    @views (dataset.domain[1:(end-1)] .+ dataset.domain[2:end]) ./ 2

make_model_domain(::ContiguouslyBinned, dataset::BinnedData) = dataset.domain
make_objective(::ContiguouslyBinned, dataset::BinnedData) = dataset.codomain

make_model_domain(::OneToOne, dataset::BinnedData) = dataset.domain[1:(end-1)]
make_objective(::OneToOne, dataset::BinnedData) = dataset.codomain

make_output_domain(layout::AbstractLayout, dataset::BinnedData) =
    make_model_domain(layout, dataset)

function make_objective_variance(::AbstractLayout, dataset::BinnedData{V})::V where {V}
    if !isnothing(dataset.codomain_variance)
        dataset.codomain_variance
    else
        # TODO: i dunno just something
        ones(eltype(V), length(dataset.codomain))
    end
end

function objective_transformer(::AbstractLayout, dataset::BinnedData)
    function _transformer!!(domain, objective)
        @views objective
    end
    function _transformer!!(output, domain, objective)
        @. output = objective
    end
    _transformer!!
end

make_label(dataset::BinnedData) = dataset.name

function _printinfo(io::IO, data::BinnedData)
    println(
        io,
        Crayons.Crayon(foreground = :cyan),
        "BinnedData",
        Crayons.Crayon(reset = true),
        " with ",
        Crayons.Crayon(foreground = :cyan),
        length(data.domain),
        Crayons.Crayon(reset = true),
        " data points:",
    )
    dmin, dmax = prettyfloat.(extrema(data.domain))
    comin, comax = prettyfloat.(extrema(data.codomain))
    descr = """  Name                  : $(data.name)
      . Domain (min/max)    : ($dmin, $dmax)
      . Codomain (min/max)  : ($comin, $comax)
    """
    print(io, descr)
end

export BinnedData
