
@with_kw struct PowerLaw <: AbstractSpectralModel{Additive}
    K::Float64 = 1.0
    a::Float64 = 0.5
end

# this is what i would like the implementation to be
# function invokemodel!(output, err, energy, m::PowerLaw)
#     @. output = energy ^ (-m.a)
# end

# this is what it has to be right now
function invokemodel!(output, err, energy, m::PowerLaw)
    alpha = 1 - m.a
    alpha_inv = inv(alpha)
    @. output = alpha_inv * energy^alpha
    # this is the normalisation needed to ensure flux is
    # integrable over energy bins
    @views output[1:end-1] .= output[2:end] .- output[1:end-1]
    output
end

export PowerLaw
