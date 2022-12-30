using Test
using SpectralFitting

include("../dummies.jl")

model = DummyAdditive()

# check parameter conversion
free_params = [2.0, 2.0]
res = SpectralFitting.free_parameters_to_named_tuple(free_params, model)
@test res == (; K = free_params[1], a = free_params[2])

# check rebuild
outmodel = SpectralFitting.remake_with_number_type(model)
@test outmodel.K == 1.0
@test outmodel.a == 1.0
@test outmodel.b == 5.0

# check free parameter conversion
outmodel = SpectralFitting.remake_with_free(model, free_params)
@test outmodel.K == free_params[1]
@test outmodel.a == free_params[2]
@test outmodel.b == 5.0

# check destructuring
out = SpectralFitting.all_parameters_to_named_tuple(
    SpectralFitting.remake_with_number_type(model),
)
@test out == (; K = 1.0, a = 1.0, b = 5.0)

# check it works with closures okay too
model = DummyAdditiveTableModel()
res = SpectralFitting.free_parameters_to_named_tuple(free_params, model)
@test res == (; K = free_params[1], a = free_params[2])

outmodel = SpectralFitting.remake_with_number_type(model)
@test outmodel.K == 1.0
@test outmodel.a == 1.0
@test outmodel.b == 2.0
