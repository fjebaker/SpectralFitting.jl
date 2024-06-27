function _convolve_1d_same_domain!(
    output::Vector{T},
    vec_A::Vector{T},
    kernel::Vector{T},
) where {T<:Real}
    # Based on https://discourse.julialang.org/t/97658/15
    @assert length(output) == length(vec_A)
    @assert length(output) == length(kernel)

    fill!(output, 0)

    @turbo for i in eachindex(output)
        total = zero(T)
        for k in eachindex(output)
            ib0 = (i >= k)
            oa = ib0 ? vec_A[i-k+1] : zero(T)
            total += kernel[k] * oa
        end
        output[i] = total
    end
    output
end

convolve!(output, A, kernel) = _convolve_1d_same_domain!(output, A, kernel)
function convolve(A, kernel)
    output = similar(A)
    convolve!(output, A, kernel)
end
