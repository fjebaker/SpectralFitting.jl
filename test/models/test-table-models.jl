using Test
using SpectralFitting

include("../dummies.jl")


# multiplicative

model = DummyMultiplicativeTableModel()

# ensure the symbols are being parsed correctly
all_symbols = SpectralFitting.all_parameter_symbols(model)
@test all_symbols == (:a, :b)
free_symbols = SpectralFitting.free_parameter_symbols(model)
@test free_symbols == (:a,)
frozen_symbols = SpectralFitting.frozen_parameter_symbols(model)
@test frozen_symbols == (:b,)


# can we invoke the table models alright
energy = collect(range(0.1, 100.0, 100))
out_flux = invokemodel(energy, model)
@test all(out_flux .== 2)

# and compose them
cm = model * model
out_flux = invokemodel(energy, cm)
@test all(out_flux .== 4.0)

# additive

model = DummyAdditiveTableModel()

all_symbols = SpectralFitting.all_parameter_symbols(model)
@test all_symbols == (:K, :a, :b)
free_symbols = SpectralFitting.free_parameter_symbols(model)
@test free_symbols == (:K, :a)
frozen_symbols = SpectralFitting.frozen_parameter_symbols(model)
@test frozen_symbols == (:b,)

out_flux = invokemodel(energy, model)
@test all(out_flux .== 3.0)

# composition
cm = model + model
out_flux = invokemodel(energy, cm)
@test all(out_flux .== 6.0)

# compose multiplicative with additive
cm = DummyMultiplicativeTableModel() * DummyAdditiveTableModel()
out_flux = invokemodel(energy, cm)
@test all(out_flux .== 6.0)