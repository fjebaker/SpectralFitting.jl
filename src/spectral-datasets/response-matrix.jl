export ResponseMatrix, fold_response

function ResponseMatrix(fits, ::AbstractMission, ::Type{T})::ResponseMatrix{T} where {T}
    # temporary structs to make parsing easier
    rmf = OGIP_RMF_Matrix(fits, T)
    chan = OGIP_RMF_Channels(fits, T)
    # allocate sparse matrix
    matrix = spzeros(T, rmf.number_of_channels, rmf.number_of_energies)
    # populate it
    build_matrix_response!(matrix, rmf)
    # unpack and return
    ResponseMatrix(
        matrix,
        chan.channels,
        chan.energy_bins_low,
        chan.energy_bins_high,
        rmf.energy_bins_low,
        rmf.energy_bins_high,
    )
end

Base.:*(rm::ResponseMatrix, flux) = fold_response(flux, rm)

@fastmath fold_response(flux, rm::ResponseMatrix) = rm.matrix * flux

function fold_response(flux, energy, rm::ResponseMatrix)
    if length(rm.energy_bins_low) != length(energy) - 1
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
    energy_low = zeros(T, N)
    energy_high = zeros(T, N)

    grouping_indices_callback(indices) do (i, index1, index2)
        @views R[i, :] .= vec(sum(rm.matrix[index1:index2, :], dims = 1))
        energy_low[i] = rm.channel_energy_bins_low[index1]
        energy_high[i] = rm.channel_energy_bins_high[index1]
    end

    ResponseMatrix(
        R,
        collect(1:N),
        energy_low,
        energy_high,
        rm.energy_bins_low,
        rm.energy_bins_high,
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
