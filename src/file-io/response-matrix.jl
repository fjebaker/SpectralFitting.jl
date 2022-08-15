export ResponseMatrix, fold_response, energy_bins

struct ResponseMatrix{T}
    matrix::SparseMatrixCSC{T,Int}
    channels::Vector{Int}
    low_energy_bins::Vector{T}
    high_energy_bins::Vector{T}
end

# folding response

@fastmath fold_response(flux, rm::ResponseMatrix) = rm.matrix * flux

function fold_response(flux, energy, rm::ResponseMatrix)
    if length(rm.low_energy_bins) != length(energy) - 1
        out_flux = rebin_flux(flux, energy, rm)
        fold_response(out_flux, rm)
    else
        fold_response(flux, rm)
    end
end

get_energy_bins(rm::ResponseMatrix{T}) where {T} = get_energy_bins(rm, T)

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, rm::ResponseMatrix{T}) where {T}
    nchans = length(rm.channels)
    println(io, "ResponseMatrix with $nchans channels:")
    Base.print_array(io, rm.matrix)
end

function group_response!(rm::ResponseMatrix, indices, T)
    N = length(indices) - 1
    R = spzeros(T, (N,N))
    energy_low = zeros(T, N)
    energy_high = zeros(T, N)
    channels = collect(1:N)

    grouping_indices_callback(indices) do (i, index1, index2)
        energy_low[i] = rm.low_energy_bins[index1]
        energy_high[i] = rm.high_energy_bins[index2]

        R[i, :] = @views sum(rm.matrix[index1:index2, 1:N]; dims=1)
    end

    ResponseMatrix(R, channels, energy_low, energy_high)
end
