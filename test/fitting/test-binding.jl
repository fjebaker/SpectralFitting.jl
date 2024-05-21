using SpectralFitting, Test

include("../dummies.jl")

dummy_data1 = make_dummy_dataset((E) -> E^(-3.0) + E^(-2.3); units = u"counts / (s * keV)")

model1 = PowerLaw() + PowerLaw()
model2 = sum(PowerLaw() for i = 1:3)

prob = FittingProblem(model1 => dummy_data1, model2 => dummy_data1)

_, mapping = SpectralFitting._build_parameter_mapping(prob.model, prob.bindings)
@test mapping == ([1, 2, 3, 4], [5, 6, 7, 8, 9, 10])

# now we bind a parameter
bind!(prob, :K_1)

_, mapping = SpectralFitting._build_parameter_mapping(prob.model, prob.bindings)
@test mapping == ([1, 2, 3, 4], [1, 5, 6, 7, 8, 9])


model1 = PowerLaw()
model2 = PowerLaw() + PowerLaw()

# binding weird model combinations
prob = FittingProblem(model1 => dummy_data1, model2 => dummy_data1)
bind!(prob, 1 => :K, 2 => :K_2)

_, mapping = SpectralFitting._build_parameter_mapping(prob.model, prob.bindings)
@test mapping == ([1, 2], [3, 4, 1, 5])

model1 =
    PowerLaw(K = FitParam(3.0), a = FitParam(4.0)) +  # K_2 a_2   3 4
    PowerLaw(K = FitParam(1.0), a = FitParam(2.0))    # K_1 a_1   1 2 
model2 =
    PowerLaw(K = FitParam(9.0), a = FitParam(10.0)) + # K_3 a_3   5 6
    PowerLaw(K = FitParam(7.0), a = FitParam(8.0)) +  # K_2 a_2   3 4 
    PowerLaw(K = FitParam(5.0), a = FitParam(6.0))    # K_1 a_1   1 2

# parameter array that is constructed -> [ 1 2 3 4  1 2 3 4 5 6]
@test keys(SpectralFitting.parameter_named_tuple(model1)) == (:K_1, :a_1, :K_2, :a_2)
@test keys(SpectralFitting.parameter_named_tuple(model2)) ==
      (:K_1, :a_1, :K_2, :a_2, :K_3, :a_3)

prob = FittingProblem(model1 => dummy_data1, model2 => dummy_data1)
bind!(prob, :K_1, :K_2)
@test prob.bindings == [[1 => 1, 2 => 1], [1 => 3, 2 => 3]]

parameters, mapping = SpectralFitting._build_parameter_mapping(prob.model, prob.bindings)
@test mapping == ([1, 2, 3, 4], [1, 5, 3, 6, 7, 8])

values = trunc.(Int, get_value.(parameters))
# should have removed K_1 and K_2 from model 2
@test values == [1, 2, 3, 4, 6, 8, 9, 10]

model1 = PowerLaw()
model2 = PowerLaw() + PowerLaw()

# test multiple bindings with multiple models and multiple datasets
prob = FittingProblem(model1 => dummy_data1, model1 => dummy_data1, model1 => dummy_data1)
bind!(prob, :K, :a)
_, mapping = SpectralFitting._build_parameter_mapping(prob.model, prob.bindings)
@test mapping == ([1, 2], [1, 2], [1, 2])

# 3 models, 1 bound parameter
prob = FittingProblem(model1 => dummy_data1, model1 => dummy_data1, model1 => dummy_data1)
bind!(prob, :a)
_, mapping = SpectralFitting._build_parameter_mapping(prob.model, prob.bindings)
@test mapping == ([1, 2], [3, 2], [4, 2])

# 3 models (2 the same, 1 different), bind one parameter between the two
prob = FittingProblem(model1 => dummy_data1, model2 => dummy_data1, model1 => dummy_data1)
bind!(prob, 1 => :a, 2 => :a_2, 3 => :a)
_, mapping = SpectralFitting._build_parameter_mapping(prob.model, prob.bindings)
@test mapping == ([1, 2], [3, 4, 5, 2], [6, 2])

# 2 models, both different, bind a parameter that is only in one model (check it does no-ops okay)
# prob = FittingProblem(model1 => dummy_data1, model2 => dummy_data1)
# bind!(prob, :K)
# note that this does not work at present because `_get_index_of_symbol` throws an error if the symbol is not found
