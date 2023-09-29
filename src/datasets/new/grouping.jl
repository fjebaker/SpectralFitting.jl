struct GroupingIterator
    grouping::Vector{Int}
end

function Base.length(grouping::GroupingIterator)
    count(==(1), grouping.grouping)
end

const GroupingIteratorItem = Union{Nothing,Tuple{Int,Int,Int}}
const GroupingIteratorState = Tuple{Int,Union{Nothing,Int}}
const GroupingIteratorItemAndState = Tuple{GroupingIteratorItem,GroupingIteratorState}

function Base.iterate(grouping::GroupingIterator)
    i1 = findfirst(==(1), grouping.grouping)
    if isnothing(i1)
        return nothing
    end
    iterate(grouping, (0, i1))
end

function Base.iterate(grouping::GroupingIterator, state)
    i::Int, i1::Int = state
    if (i >= lastindex(grouping.grouping)) || (i1 >= lastindex(grouping.grouping))
        return nothing
    end
    index = findnext(==(1), grouping.grouping, i1 + 1)
    i2::Int = if isnothing(index)
        lastindex(grouping.grouping)
    else
        index - 1
    end
    (i + 1, i1, i2), (i + 1, i2 + 1)
end

function regroup!(vector::Vector, grouping::Vector{Int})
    itt = GroupingIterator(grouping)
    for grp in itt
        regroup_vector!(vector, grp)
    end
    resize!(vector, length(itt))
    vector
end

regroup_vector!(vector::Vector, grouping::Tuple{Int,Int,Int}) =
    vector[grouping[1]] = @views sum(vector[grouping[2]:grouping[3]])
regroup_vector!(output::Vector, vector::Vector, grouping::Tuple{Int,Int,Int}) =
    output[grouping[1]] = @views sum(vector[grouping[2]:grouping[3]])

# function _regroup_warnings(data::SpectralDataset)
#     if !all(==(0), data.spectrum.quality[data.mask])
#         @warn "Dataset contains bad channels, and regrouping may be erroneous. Proceed only if you know what you're doing.\nBad channels masked out by using `mask_bad_channels!`."
#     end
# end

const BAD_QUALITY = 1
const GOOD_QUALITY = 0
function regroup_quality_vector!(quality::Vector{Int}, grouping::Tuple{Int,Int,Int})
    qs = @views quality[grouping[2]:grouping[3]]
    if any(==(BAD_QUALITY), qs)
        quality[grouping[1]] = BAD_QUALITY
    else
        quality[grouping[1]] = GOOD_QUALITY
    end
end

export regroup!


### old

# function regroup!(data::SpectralDataset{U,M,T}, grouping) where {U,M,T}
#     _regroup_warnings(data)

#     indices = grouping_to_indices(grouping)
#     N = length(indices) - 1

#     energy_low = zeros(T, N)
#     energy_high = zeros(T, N)
#     new_mask = zeros(Bool, N)

#     old_bins_low = data.bins_low
#     old_bins_high = data.bins_high

#     _grp1, _fin1 = regroup_callbacks(data.spectrum, N, data.units)
#     _grp2, _fin2 = regroup_callbacks(data.background, N, data.units)

#     grouping_indices_callback(indices) do (i, index1, index2)
#         energy_low[i] = old_bins_low[index1]
#         energy_high[i] = old_bins_high[index2]
#         # if not all are masked out, set true
#         new_mask[i] = any(==(true), @views(data.mask[index1:index2]))
#         _grp1(i, index1, index2)
#         _grp2(i, index1, index2)
#     end

#     response = group_response_channels(data.response, grouping)

#     data.bins_low = energy_low
#     data.bins_high = energy_high
#     data.spectrum = _fin1(indices)
#     data.background = _fin2(indices)
#     data.response = response
#     data.mask = new_mask
#     data
# end

# function regroup_callbacks(spectrum::Spectrum{T}, N, units) where {T}
#     channels_grouped = zeros(Int, N)
#     values_grouped = zeros(T, N)

#     errors_grouped = ismissing(spectrum.errors) ? missing : zeros(T, N)
#     function _grouper(i, index1, index2)
#         channels_grouped[i] = spectrum.channels[index1]

#         vs = @views sum(spectrum.values[index1:index2])
#         values_grouped[i] = vs

#         # TODO: rework errors
#         if !ismissing(errors_grouped)
#             # errors are calculated on **counts** not on the rates
#             if units == SpectralUnits.u"counts"
#                 errors_grouped[i] = count_error(vs, 1.0)
#             elseif units == SpectralUnits.u"counts / s"
#                 vs = vs * spectrum.exposure_time
#                 es = count_error(vs, 1.0)
#                 errors_grouped[i] = es / spectrum.exposure_time
#             else
#                 error("No method for grouping errors with given spectral units ($units).")
#             end
#         end
#     end
#     function _finalize(indices)
#         quality_grouped = group_quality_vector(spectrum.quality, indices)
#         grouping = ones(Int, N)
#         Spectrum(
#             channels_grouped,
#             quality_grouped,
#             grouping,
#             values_grouped,
#             spectrum.unit_string,
#             spectrum.exposure_time,
#             spectrum.background_scale,
#             spectrum.area_scale,
#             spectrum.error_statistics,
#             errors_grouped,
#             spectrum.systematic_error,
#             spectrum.telescope,
#             spectrum.instrument,
#         )
#     end
#     _grouper, _finalize
# end
