export SimpleDataset

struct SimpleDataset{VecType,XUnitType,YUnitType,ErrType} <: AbstractDataset
    x::VecType
    y::VecType
    x_units::XUnitType
    y_units::YUnitType
    y_err::ErrType
    x_err::ErrType
    name::String
    function SimpleDataset(
        name::String,
        x::V,
        y::V;
        xerr = nothing,
        yerr = nothing,
        y_units::YU = SpectralUnits.u"counts",
        x_units::XU = SpectralUnits.u"eV",
    ) where {V,YU,XU}
        new{V,XU,YU,typeof(xerr)}(x, y, x_units, y_units, xerr, yerr, name)
    end
end

target_vector(data::SimpleDataset) = data.y
domain_vector(data::SimpleDataset) = data.x
# if no standard deviations given
variance_vector(data::SimpleDataset{V,X,Y,<:Nothing}) where {V,X,Y} =
    ones(eltype(data.y), size(data.y))
# else standard procedure
variance_vector(data::SimpleDataset) = data.yerr .^ 2
