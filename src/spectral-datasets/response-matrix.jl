export ResponseMatrix, fold_response

function ResponseMatrix(rmf::OGIP_RMF{T}) where {T}
    rm = rmf.ogip_matrix
    chan = rmf.ogip_rmf_channels
    # allocate sparse matrix
    matrix = spzeros(T, rm.number_of_channels, rm.number_of_energies)
    # populate it
    build_matrix_response!(matrix, rm)
    # unpack and return 
    ResponseMatrix(
        matrix,
        chan.channels,
        chan.bins_low,
        chan.bins_high,
        rm.bins_low,
        rm.bins_high,
    )
end

Base.:*(rm::ResponseMatrix, flux) = fold_response(flux, rm)

@fastmath fold_response(flux, rm::ResponseMatrix) = rm.matrix * flux

function fold_response(flux, energy, rm::ResponseMatrix)
    if length(rm.bins_low) != length(energy) - 1
        out_flux = rebin_flux(flux, energy, rm)
        fold_response(out_flux, rm)
    else
        fold_response(flux, rm)
    end
end

function group_response_channels(rm::ResponseMatrix{T}, grouping) where {T}
    indices = grouping_to_indices(grouping)
    N = length(indices) - 1
    R = spzeros(T, (N, size(rm.matrix, 2)))
    bin_low = zeros(T, N)
    bin_high = zeros(T, N)

    grouping_indices_callback(indices) do (i, index1, index2)
        @views R[i, :] .= vec(sum(rm.matrix[index1:index2, :], dims = 1))
        bin_low[i] = rm.channel_bins_low[index1]
        bin_high[i] = rm.channel_bins_high[index1]
    end

    ResponseMatrix(
        R,
        collect(1:N),
        bin_low,
        bin_high,
        rm.bins_low,
        rm.bins_high,
    )
end

function build_matrix_response!(R, rmf::OGIP_RMF_Matrix)
    for (i, (F, N)) in enumerate(zip(eachcol(rmf.Fchan), eachcol(rmf.Nchan)))
        M = rmf.matrix_rows[i]
        index = 1
        for (first, len) in zip(F, N)
            if len == 0
                break
            end
            first -= rmf.first_channel
            @views R[first+1:first+len, i] .= M[index:index+len-1]
            index += len
        end
    end
end

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

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, rm::ResponseMatrix{T}) where {T}
    nchans = length(rm.channels)
    println(io, "ResponseMatrix with $nchans channels:")
    Base.print_array(io, rm.matrix)
end
