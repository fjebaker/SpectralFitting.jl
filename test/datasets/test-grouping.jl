using SpectralFitting, Test

grouping = [1, 0, 0, 0, 1, 0, 0, 1, 1, 0]

items = [i for i in SpectralFitting.GroupingIterator(grouping)]
@test items == [(1, 1, 4), (2, 5, 7), (3, 8, 8), (4, 9, 10)]

data = collect(range(0.0, 5.0, 10))
SpectralFitting.regroup!(data, grouping)

# should modify inplace
@test data ≈ [3.3333333333333335, 8.333333333333334, 3.888888888888889, 9.444444444444445]
