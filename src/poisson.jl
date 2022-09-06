# poisson statistics and related utility functionso

"""
    count_error(k, σ)

Gives the error on `k` mean photon counts to a significance of `σ` standard deviations 
about the mean.

Derived from likelihood of binomial distributions being the beta function.
"""
function count_error(k, σ)
    p = cdf(Normal(), σ)
    kₑ = gamma_inc_inv(k+1, p, 1-p)
    abs(k - kₑ)
end