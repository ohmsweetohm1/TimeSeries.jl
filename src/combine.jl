import Base: merge, vcat, map

###### merge ####################

function merge{T,N,M,D}(ta1::TimeArray{T,N,D}, ta2::TimeArray{T,M,D},
                              method::Symbol=:inner; colnames::Vector=[], meta::Any=Void)

    ta1.meta == ta2.meta && meta == Void && (meta = ta1.meta)

    if method == :inner

        idx1, idx2 = overlaps(ta1.timestamp, ta2.timestamp)
        vals = [ta1[idx1].values ta2[idx2].values]
        ta = TimeArray(ta1[idx1].timestamp, vals, [ta1.colnames; ta2.colnames], meta)

    elseif method == :left

        new_idx2, old_idx2 = overlaps(ta1.timestamp, ta2.timestamp)
        right_vals = NaN * zeros(length(ta1), length(ta2.colnames))
        right_vals[new_idx2, :]  = ta2.values[old_idx2, :]
        ta = TimeArray(ta1.timestamp, [ta1.values right_vals], [ta1.colnames; ta2.colnames], meta)

    elseif method == :right

        ta = merge(ta2, ta1, :left)
        ncol2 = length(ta2.colnames)
        vals = [ta.values[:, (ncol2+1):end] ta.values[:, 1:ncol2]]
        ta = TimeArray(ta.timestamp, vals, [ta1.colnames; ta2.colnames], meta)

    elseif method == :outer

        timestamps = sorted_unique_merge(ta1.timestamp, ta2.timestamp)
        ta = TimeArray(timestamps, zeros(length(timestamps), 0), UTF8String[], Void)
        ta = merge(ta, ta1, :left)
        ta = merge(ta, ta2, :left, meta=meta)

    else
        error("merge method must be one of :inner, :left, :right, :outer")
    end

    return setcolnames!(ta, colnames)

end

# collapse ######################

function collapse{T,N,D}(ta::TimeArray{T,N,D}, f::Function; period::Function=week)

    w = [period(ta.timestamp[t]) for t in 1:length(ta)] # get weekly id from entire array
    z = Int[]; j = 1
    for i=1:length(ta) - 1 # create unique period ID array
        if w[i] != w[i+1]
            push!(z, j)
            j = j+1
        else
            push!(z,j)
        end         
    end
    
    # account for last row
    w[length(ta)]  ==  w[length(ta)-1] ? # is the last row the same period as 2nd to last row?
    push!(z, z[length(z)]) :  
    push!(z, z[length(z)] + 1)  
   
    # pre-allocate timestamp and value arrays
    tstamps = Array{D}(maximum(z))
    vals    = zeros(maximum(z)) # number of unique periods
    #replace their values except for the last row 
    for i = 1:maximum(z)-1  # iterate over period ID groupings
        temp       = ta[findfirst(z .== i):findfirst(z .== i+1)-1] # period's worth that will be squished
        tstamps[i] = temp[length(temp)].timestamp[1]
        vals[i]    = f(temp.values)
    end
    # and once again account for the last temp that isn't looped on above
    lasttemp = ta[findfirst(z .== maximum(z)):length(ta)]
    lasttempindex = lasttemp[length(lasttemp)].timestamp
    lasttempvalue = f(ta.values)
    
    # complete the tstamp and vals arrays with the last value
    tstamps[length(tstamps)] = lasttempindex[1]
    vals[length(vals)]       = lasttempvalue

    TimeArray(tstamps, vals, ta.colnames, ta.meta)
end

# vcat ######################

function vcat{T,N,D}(TA::TimeArray{T,N,D}...)
    # Check all meta fields are identical. 
    prev_meta = TA[1].meta
    for ta in TA
        if ta.meta != prev_meta
            error("metadata doesn't match")
        end
    end
    
    # Check column names are identical. 
    prev_colnames = TA[1].colnames
    for ta in TA
        if ta.colnames != prev_colnames
            error("column names don't match")
        end
    end
    
    # Concatenate the contents. 
    timestamps = vcat([ta.timestamp for ta in TA]...)
    values = vcat([ta.values for ta in TA]...)
    
    order = sortperm(timestamps)
    if length(TA[1].colnames) == 1 # Check for 1D to ensure values remains a 1D vector. 
        return TimeArray(timestamps[order], values[order], TA[1].colnames, TA[1].meta)
    else
        return TimeArray(timestamps[order], values[order, :], TA[1].colnames, TA[1].meta)
    end
end

# map ######################

function map{T,N,D,A}(f::Function, ta::TimeArray{T,N,D,A})
    timestamps = similar(ta.timestamp)
    values = similar(ta.values)
    
    for i in 1:length(ta)
        timestamps[i], values[i, :] = f(ta.timestamp[i], vec(ta.values[i, :]))
    end
    
    order = sortperm(timestamps)
    if length(ta.colnames) == 1 # Check for 1D to ensure values remains a 1D vector. 
        return TimeArray(timestamps[order], values[order], ta.colnames, ta.meta)
    else
        return TimeArray(timestamps[order], values[order, :], ta.colnames, ta.meta)
    end
end
