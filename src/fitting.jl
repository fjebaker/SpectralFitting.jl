export fitparams!

function fitparams!(
    psm::ProcessedSpectralModel,
    rmf::ResponseMatrix,
    grp,
    energy = energybins(rmf);
    kwargs...,
)
    target = grp.RATE ./ grp.E_DELTA
    error_target = grp.STAT_ERR ./ grp.E_DELTA
    __fitparams!(psm, rmf, target, error_target, energy, grp.CHANNEL; kwargs...)
end

function __fitparams!(
    psm::ProcessedSpectralModel,
    rmf::ResponseMatrix,
    target,
    error_target,
    energy,
    channels;
    kwargs...,
)
    model, params = build_simple(psm)
    fit = __lsq_fit(model, params, rmf, target, error_target, energy, channels; kwargs...)

    means = LsqFit.coef(fit)
    errs = LsqFit.stderror(fit)

    for (p, m, e) in zip(filter(!isfrozen, last.(psm.parameters)), means, errs)
        p.val = m
        p.lower_bound = m - 2e
        p.upper_bound = m + 2e
    end

    fit
end

function __lsq_fit(model, params, rmf, target, error_target, energy, channels; kwargs...)
    p0 = value.(params)

    @show p0
    @show size(energy)

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
