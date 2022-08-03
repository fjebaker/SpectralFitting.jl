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
    fit = __lsq_fit(model, get_value.(params), get_lowerlimit.(params), get_upperlimit.(params), rmf, target, error_target, energy, channels; kwargs...)

    means = LsqFit.coef(fit)
    errs = LsqFit.stderror(fit)

    for (p, m, e) in zip(params, means, errs)
        set_value!(p, m)
        set_error!(p, e)
    end

    fit
end

@noinline function __lsq_fit(model, p0, lb, ub, rmf, target, error_target, energy, channels; kwargs...)
    fluxes = makefluxes(energy)
    foldedmodel(x, p) = @views foldresponse(rmf, model(fluxes..., x, p), x)[channels]

    # seems to cause nans and infs
    cov_err = 1 ./ (error_target .^ 2)

    fit = LsqFit.curve_fit(
        foldedmodel,
        energy,
        target,
        cov_err,
        p0;
        lower = lb,
        upper = ub,
        kwargs...,
    )
    fit
end
