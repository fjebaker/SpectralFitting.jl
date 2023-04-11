function group_quality_vector(quality, inds)
    new_quality = zeros(Int, length(inds) - 1)
    # bad_channels = Int[]
    grouping_indices_callback(inds) do (i, index1, index2)
        qs = @views quality[index1:index2]
        # ms = @views mask[index1:index2]
        if any(!=(0), qs) # if any are bad quality
            new_quality[i] = 1
            # bad channel if not all of them are masked out
            # if !all(==(false), ms)
            #     push!(bad_channels, i)
            # end
        else
            new_quality[i] = 0
        end
    end
    # if !isempty(bad_channels)
    #     warn_bad_channels(bad_channels)
    # end
    new_quality
end

function warn_bad_channels(bad_channels)
    warns = String[]
    N = length(bad_channels)
    j = 1
    for (i, Δ) in enumerate(diff(bad_channels))
        if Δ != 1
            if i - j == 0
                push!(warns, "$i")
            else
                push!(warns, "$j-$i")
            end
            j = i
        end
    end
    if N - j == 0
        push!(warns, "$j")
    else
        push!(warns, "$j-$N")
    end
    chans = join(warns, ", ")
    @warn "Grouped channels $chans contain bad quality channels."
end

infer_units(s) = SpectralUnits.infer_units(s)

# include missions
include("no-associated-mission.jl")
include("xmm-newton-mission.jl")
