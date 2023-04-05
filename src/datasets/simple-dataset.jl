export SimpleDataset

struct SimpleDataset{VecType,XUnitType,YUnitType,ErrType} <: AbstractDataset
    x::VecType
    y::VecType
    x_units::XUnitType
    y_units::YUnitType
    x_err::ErrType
    y_err::ErrType
    name::String
    function SimpleDataset(
        name::String,
        x::V,
        y::V;
        x_err = nothing,
        y_err = nothing,
        y_units::YU = SpectralUnits.u"counts",
        x_units::XU = SpectralUnits.u"eV",
    ) where {V,YU,XU}
        _x_err = isnothing(x_err) ? ones(eltype(x), size(x)) : x_err
        _y_err = isnothing(y_err) ? ones(eltype(y), size(y)) : y_err

        # check error dimensions
        if size(_x_err) != size(x)
            error("Size of x_err must be the same as x ($(sizeof(_x_err)) != $(sizeof(x)))")
        end
        if size(_y_err) != size(y)
            error("Size of y_err must be the same as y ($(sizeof(_y_err)) != $(sizeof(y)))")
        end

        new{V,XU,YU,typeof(_x_err)}(x, y, x_units, y_units, _x_err, _y_err, name)
    end
end

target_vector(data::SimpleDataset) = data.y
domain_vector(data::SimpleDataset) = data.x
target_variance(data::SimpleDataset) = data.y_err .^ 2
