export fitparams!

function fitparams!(
    psm::ProcessedSpectralModel,
    dataset::AbstractSpectralDataset,
    energy = rmfenergybins(dataset);
    kwargs...,
)
    model, params = build_simple(psm)
    fitparams!(params, model, dataset, energy; kwargs...)
end

function fitparams!(
    params,
    model,
    dataset::AbstractSpectralDataset,
    energy = rmfenergybins(dataset);
    kwargs...,
)
    __fitparams!(
        params,
        model,
        response(dataset),
        countbins(dataset),
        counterrors(dataset),
        energy,
        channels(dataset);
        kwargs...,
    )
end

function __fitparams!(
    params,
    model,
    rmf::ResponseMatrix,
    target,
    error_target,
    energy,
    channels;
    kwargs...,
)
    fit = __lsq_fit(model, params, rmf, target, error_target, energy, channels; kwargs...)

    means = LsqFit.coef(fit)
    errs = LsqFit.stderror(fit)

    for (p, m, e) in zip(params, means, errs)
        setvalue!(p, m)
        setlowerbound!(p, m - 2e)
        setupperbound!(p, m + 2e)
    end

    fit
end

function __lsq_fit(model, params, rmf, target, error_target, energy, channels; kwargs...)
    p0 = value.(params)

    fluxes = makefluxes(energy)
    foldedmodel(x, p) = @views foldresponse(rmf, model(fluxes..., x, p), x)[channels]

    # seems to cause nans and infs
    # cov_err = 1 ./ (error_target .^ 2)

    fit = LsqFit.curve_fit(
        foldedmodel,
        energy,
        target,
        # cov_err,
        p0;
        lower = lowerbound.(params),
        upper = upperbound.(params),
        kwargs...,
    )
    fit
end
