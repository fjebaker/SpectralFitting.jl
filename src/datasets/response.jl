# TODO: could be Response or Redistribution : how do we track this? 
mutable struct ResponseMatrix{T}
    matrix::SparseMatrixCSC{T,Int}
    channels::Vector{Int}
    channel_bins_low::Vector{T}
    channel_bins_high::Vector{T}
    bins_low::Vector{T}
    bins_high::Vector{T}
end

function regroup!(resp::ResponseMatrix{T}, grouping) where {T}
    itt = GroupingIterator(grouping)
    new_matrix = zeros(T, size(resp.matrix))
    for grp in itt
        slice = @views resp.matrix[grp[2]:grp[3], :]
        new_matrix[grp[1], :] .= sum(slice, dims = 1) |> vec

        resp.channels[grp[1]] = resp.channels[grp[2]]

        resp.channel_bins_low[grp[1]] = resp.channel_bins_low[grp[2]]
        resp.channel_bins_high[grp[1]] = resp.channel_bins_high[grp[2]]
    end
    resp.matrix = sparse(new_matrix)
    resize!(resp, length(itt))
end

function Base.resize!(response::ResponseMatrix, n)
    response.matrix = response.matrix[1:n, :]
    resize!(response.channels, n)
    resize!(response.channel_bins_low, n)
    resize!(response.channel_bins_high, n)
end

# todo: currently unused; do we need this?
function normalise_rows!(matrix)
    # returns weights
    @views map(1:size(matrix, 2)) do
        w = sum(matrix[i, :])
        if ΣR > 0.0
            matrix[i, :] .= matrix[i, :] / w
        end
        w
    end
end

function build_response_matrix!(
    R,
    f_chan::Matrix,
    n_chan::Matrix,
    matrix_rows::Vector,
    first_channel,
)
    for (i, (F, N)) in enumerate(zip(eachcol(f_chan), eachcol(n_chan)))
        M = matrix_rows
        index = 1
        for (first, len) in zip(F, N)
            if len == 0
                break
            end
            first -= first_channel
            @views R[first+1:first+len, i] .= M[index:index+len-1]
            index += len
        end
    end
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, rm::ResponseMatrix{T}) where {T}
    nchans = length(rm.channels)
    println(io, "ResponseMatrix with $nchans channels:")
    Base.print_array(io, rm.matrix)
end

function _printinfo(io, resp::ResponseMatrix{T}) where {T}
    emin, emax = minimum(resp.bins_low), maximum(resp.bins_high)
    c_emin, c_emax = minimum(resp.channel_bins_low), maximum(resp.channel_bins_high)
    descr = """Reponse Matrix:
      . Channels            : $(length(resp.channels))
      . Channel E (min/max) : ($c_emin, $c_emax)
      . Matrix dims         : $(size(resp.matrix))
      . E (min/max)         : ($emin, $emax)
    """
    print(io, descr)
end

function drop_channels!(response::ResponseMatrix, inds)
    deleteat!(response.channels, inds)
    deleteat!(response.channel_bins_high, inds)
    deleteat!(response.channel_bins_low, inds)
    keep = collect(i for i = 1:size(response.matrix, 1) if i ∉ inds)
    response.matrix = response.matrix[keep, :]
    response
end

struct AncillaryResponse{T}
    bins_low::Vector{T}
    bins_high::Vector{T}
    effective_area::Vector{T}
end

regroup!(resp::AncillaryResponse, grouping) = nothing # no op

function _printinfo(io, resp::AncillaryResponse{T}) where {T}
    emin, emax = minimum(resp.bins_low), maximum(resp.bins_high)
    descr = """Ancillary Response:
      . Channels            : $(length(resp.effective_area))
      . E (min/max)         : ($emin, $emax)
    """
    print(io, descr)
end

fold_ancillary(response::ResponseMatrix, ancillary::AncillaryResponse) =
    ancillary.effective_area' .* response.matrix

fold_ancillary(response::ResponseMatrix, ::Nothing) = response.matrix

export ResponseMatrix