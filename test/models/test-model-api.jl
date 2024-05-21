using Test, SpectralFitting

include("../dummies.jl")

# single

model = DummyAdditive()
domain = collect(range(0.1, 10, 200))
output = @inferred allocate_model_output(model, domain)
@test size(output) == (199, 1)

invokemodel!(output, domain, model)
@test all(i -> isapprox(i, 6.0), output)

output = invokemodel(domain, model)
@test all(i -> isapprox(i, 6.0), output)

invokemodel!(output, domain, model)
allocated = @allocated invokemodel!(output, domain, model)
@test allocated <= 48

@test modelkind(model) == Additive()
@test numbertype(model) == Float64

# composite

model = DummyAdditive() + DummyAdditive()
output = @inferred allocate_model_output(model, domain)
@test size(output) == (199, 2)

invokemodel!(output, domain, model)
@test all(i -> isapprox(i, 12.0), output[:, 1])

allocated = @allocated invokemodel!(output, domain, model)
@test allocated <= 48

output = invokemodel(domain, model)
@test all(i -> isapprox(i, 12.0), output)

@test modelkind(model) == Additive()
@test numbertype(model) == Float64
