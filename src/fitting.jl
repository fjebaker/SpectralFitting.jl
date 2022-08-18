export fitparams!

function fitparams!(
    params,
    model,
    dataset::SpectralDataset,
    energy = energy_vector(dataset.response);
    kwargs...,
)   

    folder(flux, energy) = begin
        # fold instrument response + arf
        folded_flux = fold_response(flux, energy, dataset.response)
        # rebin and mask
        grouped_flux = @views regroup(folded_flux, dataset.meta.grouping)[dataset.mask] 
        # normalise against energy
        grouped_flux ./ dataset.energy_bin_widths
    end
    __fitparams!(
        params,
        model,
        folder,
        dataset.counts,
        dataset.countserror,
        energy,
        ;
        kwargs...,
    )
end

function __fitparams!(
    params,
    model::M,
    folder,
    target,
    error_target,
    energy
    ;
    kwargs...
) where {M<:AbstractSpectralModel}

    p0 = get_value.(params)
    lb = get_lowerlimit.(params)
    ub = get_upperlimit.(params)

    fit = __lsq_fit(
        implementation(M),
        model,
        folder,
        p0,
        lb,
        ub,
        target,
        error_target,
        energy,
        kwargs...,
    )

    means = LsqFit.coef(fit)
    errs = LsqFit.stderror(fit)

    for (p, m, e) in zip(params, means, errs)
        set_value!(p, m)
        set_error!(p, e)
    end

    fit
end

function __lsq_fit(
    ::AbstractSpectralModelImplementation,
    model,
    folder,
    p0,
    lb,
    ub,
    target,
    error_target,
    energy,
    kwargs...,
)
    frozen_p = get_value.(get_frozen_model_params(model))
    fluxes = make_fluxes(energy, flux_count(model))
    foldedmodel(x, p) = begin
        flux = invokemodel!(fluxes, x, model, p, frozen_p)
        folder(flux, x)
    end
    # seems to cause nans and infs
    cov_err = @. 1 / (error_target^2)
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

# function __lsq_fit(::JuliaImplementation, model, p0, lb, ub, rm, target, error_target, energy, mask; kwargs...)
#     frozen_p = get_value.(get_frozen_model_params(model))
#     d_fluxes = make_dual_fluxes(energy, flux_count(model))
#     foldedmodel(x, p) = begin
#         fluxes = map(i -> get_tmp(i, p[1]), d_fluxes)
#         @views fold_response(invokemodel!(fluxes, x, model, p, frozen_p), x, rm)[mask]
#     end

#     # seems to cause nans and infs
#     cov_err = 1 ./ (error_target .^ 2)

#     fit = LsqFit.curve_fit(
#         foldedmodel,
#         energy,
#         target,
#         cov_err,
#         p0;
#         lower = lb,
#         upper = ub,
#         autodiff = :forwarddiff,
#         kwargs...,
#     )
#     fit
# end
