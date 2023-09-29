using SpectralFitting, Test

spectrum = Spectrum(
    collect(1:10),
    zeros(Int, 10),
    ones(Int, 10),
    collect(range(1.0, 0.0, 10)),
    "counts",
    2.0,
    1.0,
    1.0,
    SpectralFitting.ErrorStatistics.Poisson,
    collect(range(1.0, 0.0, 10)) ./ 100,
    0.0,
    "test-telescope",
    "test-name",
)

matrix = SpectralFitting.sparse(Matrix(1.0 * SpectralFitting.I, 10, 10))
energies = collect(range(0, 15.0, 11))
response = ResponseMatrix(
    matrix,
    spectrum.channels,
    energies[1:end-1],
    energies[2:end],
    energies[1:end-1],
    energies[2:end],
)

data = SpectralData(spectrum, response)

@test make_objective(ContiguouslyBinned(), data) == [
    1.0,
    0.8888888888888888,
    0.7777777777777778,
    0.6666666666666666,
    0.5555555555555556,
    0.4444444444444444,
    0.3333333333333333,
    0.2222222222222222,
    0.1111111111111111,
    0.0,
]
@test make_domain(ContiguouslyBinned(), data) ==
      [1.5, 0.0, 1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0, 13.5]
@test make_objective_variance(ContiguouslyBinned(), data) == data.spectrum.errors .^ 2

normalize!(spectrum)

@test make_objective(ContiguouslyBinned(), data) == [
    0.5,
    0.4444444444444444,
    0.3888888888888889,
    0.3333333333333333,
    0.2777777777777778,
    0.2222222222222222,
    0.16666666666666666,
    0.1111111111111111,
    0.05555555555555555,
    0.0,
]
@test make_domain(ContiguouslyBinned(), data) ==
      [1.5, 0.0, 1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0, 13.5]
@test make_objective_variance(ContiguouslyBinned(), data) == data.spectrum.errors .^ 2
