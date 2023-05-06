using Revise
using Test
using SpectralFitting

model = PowerLaw() + PowerLaw()

SpectralFitting.FunctionGeneration.all_parameters_to_named_tuple(typeof(model))

freeparameters(model)

energy = collect(range(0.1, 10.0, 100))

invokemodel(energy, model)


SpectralFitting.FunctionGeneration.generated_model_call!(
    typeof(energy),
    typeof(energy),
    typeof(model),
    typeof(energy),
    typeof(energy),
)
