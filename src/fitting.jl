export fitparams!

function fitparams!(
    params,
    model,
    dataset::AbstractSpectralDataset,
    energy = rmfenergybins(dataset);
    kwargs...,
)
    rmf = response(dataset)
    target = countbins(dataset)
    error_target = countbins(dataset)
    channels = channels(dataset)

    frozen_p = get_frozen_model_params(model)
    fit = __lsq_fit(
        model,
        isempty(frozen_p) ? frozen_p : get_value.(frozen_p),
        get_value.(params),
        get_lowerlimit.(params),
        get_upperlimit.(params),
        rmf,
        target,
        error_target,
        energy,
        channels;
        kwargs...,
    )

    # unpack results back into parameters

    means = LsqFit.coef(fit)
    errs = LsqFit.stderror(fit)

    for (p, m, e) in zip(params, means, errs)
        set_value!(p, m)
        set_error!(p, e)
    end

    fit
end

@noinline function __lsq_fit(
    model,
    frozen_p,
    p0,
    lb,
    ub,
    rmf,
    target,
    error_target,
    energy,
    channels;
    kwargs...,
)
    fluxes = makefluxes(energy, flux_count(model))
    foldedmodel(x, p) =
        @views foldresponse(rmf, generated_model_call!(fluxes, x, model, p, frozen_p), x)[channels]

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
