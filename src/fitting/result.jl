
struct FittingResult{T,M,E,F}
    u::Vector{T}
    χ2::T
    model::M
    x::E
    folded_invoke::F
end

function Base.show(io::IO, ::MIME"text/plain", res::FittingResult)
    println(
        io,
        """FittingResult:
          Model: $(res.model)
            - Parameters  : $(res.u)
            - χ2          : $(res.χ2) 
        """,
    )
end

function unpack_lsqfit_result(res, model, f, x, y, variance)
    u = LsqFit.coef(res)
    chi2 = χ2_from_ŷyvar(f(x, u), y, variance)
    FittingResult(u, chi2, model, x, f)
end
