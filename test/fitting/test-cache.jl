using SpectralFitting, Test

# generate some fake powerlaw data
dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")

# test that both julia and xspec implementations can fit simple
model = PhotoelectricAbsorption() * PowerLaw()

config = SpectralFitting.FittingConfig(model, dummy_data)

f = SpectralFitting._f_objective(config)
params = get_value.(config.parameters)
domain = config.domain

result = f(domain, params)
# todo: currently we still allocate cus extracting the frozen parameters
# has a dynamic allocation
allocated_bytes = @allocated f(domain, params)
@test allocated_bytes < 400
