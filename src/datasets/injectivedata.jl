
struct InjectiveData{V} <: AbstractDataset
    domain::V
    codomain::V
    domain_variance::Union{Nothing,V}
    codomain_variance::Union{Nothing,V}
    name::String
    data_mask::BitVector
end

function InjectiveData(
    domain,
    codomain;
    domain_variance = nothing,
    codomain_variance = nothing,
    name = "[no-name]",
    data_mask = BitVector(fill(true, size(codomain)))
)
    InjectiveData(domain, codomain, domain_variance, codomain_variance, name, data_mask)
end

supports_contiguosly_binned(::Type{<:InjectiveData}) = true
supports_one_to_one(::Type{<:InjectiveData}) = true

function make_model_domain(::ContiguouslyBinned, dataset::InjectiveData)
    # need to expand the domain by one
    Δ = sum(diff(dataset.domain)) / (length(dataset.domain) - 1)
    domain = copy(dataset.domain)
    push!(domain, domain[end] + Δ)
    domain
end
make_objective(::ContiguouslyBinned, dataset::InjectiveData) = dataset.codomain[dataset.data_mask]

make_model_domain(::OneToOne, dataset::InjectiveData) = dataset.domain[data_mask]
make_objective(::OneToOne, dataset::InjectiveData) = dataset.codomain[data_mask]

function make_objective_variance(
    ::AbstractDataLayout,
    dataset::InjectiveData{V},
)::V where {V}
    if !isnothing(dataset.domain_variance)
        dataset.codomain_variance[dataset.data_mask]
    else
        # todo: i dunno just something
        1e-8 .* dataset.codomain[dataset.data_mask]
    end
end

objective_transformer(::AbstractDataLayout, dataset::InjectiveData) = _DEFAULT_TRANSFORMER()

make_label(dataset::InjectiveData) = dataset.name

function _printinfo(io::IO, data::InjectiveData)
    println(
        io,
        Crayons.Crayon(foreground = :cyan),
        "InjectiveData",
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

export InjectiveData
