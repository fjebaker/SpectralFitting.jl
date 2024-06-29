function _convolve_implementation!(
    output::AbstractVector{T},
    vec_A::AbstractVector{T},
    kernel::AbstractVector{T},
) where {T<:Number}
    # Based on https://discourse.julialang.org/t/97658/15
    J = length(vec_A)
    K = length(kernel)
    @assert length(output) == J + K - 1 "Ouput is $(length(output)); should be $(J + K - 1)"

    # do the kernel's side first
    for i = 1:K-1
        total = zero(T)
        for k = 1:K
            ib = (i >= k)
            oa = ib ? vec_A[i-k+1] : zero(T)
            total += kernel[k] * oa
        end
        output[i] = total
    end
    # now the middle
    for i = K:J-1
        total = zero(T)
        for k = 1:K
            oa = vec_A[i-k+1]
            total += kernel[k] * oa
        end
        output[i] = total
    end
    # and finally the end
    for i = J:(J+K-1)
        total = zero(T)
        for k = 1:K
            ib = (i < J + k)
            oa = ib ? vec_A[i-k+1] : zero(T)
            total += kernel[k] * oa
        end
        output[i] = total
    end
    output
end

convolve!(output, A, kernel) = _convolve_implementation!(output, A, kernel)
function convolve(A, kernel)
    output = zeros(eltype(A), length(A) + length(kernel) - 1)
    convolve!(output, A, kernel)
end
