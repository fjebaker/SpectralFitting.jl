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

out = SpectralFitting.all_parameters_to_named_tuple(
    SpectralFitting.remake_with_number_type(model),
)
@test out == (; K = 1.0, a = 1.0, b = 2.0)

outmodel = SpectralFitting.remake_with_number_type(model)
@test outmodel.K == 1.0
@test outmodel.a == 1.0
@test outmodel.b == 2.0


# composed models
model =
    DummyMultiplicative(a = 0.0, b = 1.0) *
    (DummyAdditive(K = 2.0, a = 3.0) + DummyAdditive(K = 4.0, a = 5.0))
out = SpectralFitting.all_parameters_to_named_tuple(model)
@test out == (;
    K_1 = 4.0,
    a_1 = 5.0,
    b_1 = 5.0,
    K_2 = 2.0,
    a_2 = 3.0,
    b_2 = 5.0,
    a_3 = 0.0,
    b_3 = 1.0,
)
