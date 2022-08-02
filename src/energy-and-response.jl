export ResponseMatrix, foldresponse, energybins, energy_to_channel, rebinflux

struct ResponseMatrix{T,E}
    matrix::SparseMatrixCSC{T,Int64}
    ebins::E
end

foldresponse(rmf::ResponseMatrix, flux) = rmf.matrix * flux
Base.:*(rmf::ResponseMatrix, flux) = rmf.matrix * flux

function energybins(rmf::ResponseMatrix{T}) where {T}
    energy = zeros(T, length(rmf.ebins.E_MAX) + 1)
    energy[1:end-1] .= rmf.ebins.E_MIN
    energy[end] = rmf.ebins.E_MAX[end]
    energy
end

function Base.show(
    io::IO,
    ::MIME{Symbol("text/plain")},
    rmf::ResponseMatrix{E,M},
) where {E,M}
    nchans = nrow(rmf.ebins)
    println(io, "ReponseMatrix with $nchans channels:")
    Base.print_array(io, rmf.matrix)
end

function energy_to_channel(energy, ebins)
    for j in eachrow(ebins)
        if (energy ≥ j.E_MIN) && (j.E_MAX > energy)
            return convert(Int64, j.CHANNEL)
        end
    end
    # for type stability
    return -1
end

rebinflux(flux, curr_energy, dest_energy::AbstractVector) =
    rebinflux(flux, curr_energy, first(dest_energy), @view(dest_energy[2:end]))
rebinflux(flux, curr_energy, ebins::DataFrame) =
    rebinflux(flux, curr_energy, first(ebins.E_MIN), ebins.E_MAX)
rebinflux(flux, curr_energy, rmf::ResponseMatrix) = rebinflux(flux, curr_energy, rmf.ebins)

function rebinflux(flux, curr_energy, first_E_min::Number, E_max_array::AbstractVector)
    N = length(E_max_array)
    out_flux = zeros(eltype(flux), N)

    # find the initial index into energy
    current_i = findfirst(>(first_E_min), curr_energy)

    for (flux_index, E_max) in enumerate(E_max_array)
        # find where energy is outside of bin
        next_i = findnext(>(E_max), curr_energy, current_i + 1)
        isnothing(next_i) && break

        # sum everything inbetween
        @views out_flux[flux_index] += sum(flux[current_i:next_i-2])

        # calculate overlap flux: width of energy bin
        ΔE = curr_energy[next_i] - curr_energy[next_i-1]
        # width of bin in current flux bin
        δE = E_max - curr_energy[next_i-1]
        ratio = (δE / ΔE)
        @views out_flux[flux_index] += ratio * flux[next_i-1]

        # if not at end of flux array, carry over
        if flux_index <= N
            out_flux[flux_index+1] += (1 - ratio) * flux[next_i-1]
        end
        current_i = next_i
    end

    out_flux
end

function foldresponse(rmf::ResponseMatrix, flux, energy)
    if nrow(rmf.ebins) != length(energy) - 1
        out_flux = rebinflux(flux, energy, rmf.ebins)
        foldresponse(rmf, out_flux)
    else
        foldresponse(rmf, flux)
    end
end
