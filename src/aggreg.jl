using Statistics, OrderedCollections, Dates, Parameters, Base.Threads

abstract type SData end

"""
InputData <: SData
Contains time series data (i.e., respecting a temporal order)
Attributes:
- region::String: optional information to specify the region data belongs to
- period_from::DateTime: date(year,month,day,hour) from which the data is from
- period_to::DateTime: date(year,month,day,hour) until which the data is
- resolution::Period: resolution of the original series (e.g., hours, days, weeks, etc.)
- nseries::Int: number of time series inputed (treated separately as if independent from each other)
- data::Dict{String,Array}: Dictionary with an entry for each attribute `[file name (attribute: e.g technology)]-[column name (node: e.g. location)]`, Each entry of the dictionary is a 2-dimensional `time-steps T x periods K`-Array holding the data
"""
mutable struct InputData
    region::String
    period_from::DateTime
    period_to::DateTime
    resolution::Period
    nseries::Int
    data::Dict{String,Array}
end

"""
SeriesInstance <: SData
Contains the series and respective attributes needed for performing the aggregation
Attributes:
- series::VecOrMat{T}: original time series
- block_size::Int: number of clusters to be aggregated at once each iteration
- stopping_k::Int: reduced number of clusters pre-established
- current_k::Int: current number of clusters
- dm::Symbol: discrepancy metric to be used to measure the difference between clusters to be aggregated
- rep_value::Symbol: function to assign a representative value for each cluster (by default, the mean of values inside the cluster)
- lseries::Int: length of series
- nseries::Int: number of series
- series_dc::VecOrMat{T}: original time series duration curves
- ord_dc::VecOrMat{S}: original order of the elements
"""
mutable struct SeriesInstance{T<:Float64,S<:Int}
    series::VecOrMat{T}
    block_size::S
    stopping_k::S
    current_k::S
    dm::Symbol
    rep_value::Symbol
    lseries::S
    nseries::S
    series_dc::VecOrMat{T}
    ord_dc::VecOrMat{S}
end

"""
ClustInstance <: SData
Contains clusters information
Fields:
- k_cent::VecOrMat{T}: 
- weights::Vector{<:Int}: 1-dimensional `periods K`-Array with the absolute weight for each period. E.g. for a yearly resolution time series of 365 days, sum(weights) = 365
- series_clust::Vector{<:Int}: clusters assignment for the original series (vector with the from-to relationship between the original series and the clusters)
- nclusters::Int: current number of clusters
- search_range::UnitRange: range in which the clusters can be merged, considering the amount of clusters and the block_size (generally equals to 1:(nclusters-block_size+1))
- dc_mode::Bool: false by default, this attribute determines if we use the duration curves order to perform the clustering
"""
@with_kw mutable struct ClustInstance{T<:Float64,S<:Int}
    k_cent::VecOrMat{T}
    weights::Vector{S}
    series_clust::Vector{S}
    nclusters::S
    search_range::UnitRange
    dc_mode::Bool = false
end

"""
AggregInstance <: SData
Contains the aggregation instance temporary values.
Attributes:
- merging_clust::UnitRange{<:Int}: 
- series_comp::VecOrMat{T}: 
- k_cent_comp::VecOrMat{T}: 
"""
mutable struct AggregInstance{T<:Float64,S<:Int}
    merging_clust::UnitRange{S}
    series_comp::VecOrMat{T}
    k_cent_comp::VecOrMat{T}
end

"""
DistUpdate <: SData
Type to store the history of minimal distances.
"""

mutable struct DistUpdate{S<:Int}
    min_dist::S
    merging_clust::UnitRange{S}
end

"""
aggreg1D(series::VecOrMat{T}, marker::Vector{Bool}, rep_value::Symbol = :mean)
The function aggregates N one-dimensional vectors concomitantly following a marker gave as the input of the function. It assumes the collection of 1D vectors is vertically concatanated so forming a 2D matrix with the first dimension equal to the size of the vectors analysed and the second dimension equal to the number of vectors analysed.
Fields:
- series::VecOrMat{T}: collection of vectors
- marker::Vector{Bool}: vector with the flag on where to aggregate - so whenever marker[i] = 1 the aggregation happens between i and i+1
- rep_value::Symbol: value to be used for representation of the aggregated interval (we use mean by default)
"""
# (first method)
function aggreg1D(series::VecOrMat{T}, marker::Vector{Bool}, representation::Symbol = :mean) where {T<:Float64}
    #### Function to aggregate vectors under a representative value, normally used in hierarchical
    #### clustering. It merges a set of vectors (series) merging the state (i) with (i + 1) whenever there's a marker

    @assert length(size(series)) <= 2 "Series to be aggregated need to be represented as a Matrix (2D)"

    ## Initialization
    lseries = size(series,1)                      # Length of vectors
    nseries = size(series,2)                      # Number of vectors
    D = 1:nseries                                 # Range for the number of series
    class = collect(1:lseries)                    # Vector to mark the new classes
    new_v = copy(series)                          # New vector to be formed after aggregation
    aux = zeros(lseries)

    ## Error control
    @assert lseries == length(marker) "Different sizes"
    @assert any((marker .!= 0) .| (marker .!= 1)) "Inconsistent marker vector"

    ##Accounts the cumulative change in clusters order
    c_change = 0                        # Change in the number of clusters
    for i in 2:lseries
        c_change -= marker[i-1]
        aux[i] = c_change
    end

    ## Update values
    class .+= aux

    ## Calculate the representative value and replace it in the new_v
    if representation == :mean
        for i in 1:class[end], d in D
            new_v[class.==i,d] .= mean(series[class.==i,d])
        end
    elseif representation == :medoid
        for i in 1:class[end], d in D
            dist = pairwise(Euclidean(),series[class.==i,d])
            s_medoid = kmedoids(dist, 1)
            s_medoid = s_medoid.medoids[1]
            new_v[class.==i,d] .= series[class.==i,d][s_medoid]
        end
    else
        @assert false "Representation method not defined."
    end

    return new_v, class
end

"""
aggreg1D(series::VecOrMat{T}, rep_value::Symbol = :mean)
The function aggregates N one-dimensional vectors concomitantly following a marker gave as the input of the function. It assumes the collection of 1D vectors is vertically concatanated so forming a 2D matrix with the first dimension equal to the size of the vectors analysed and the second dimension equal to the number of vectors analysed.
Fields:
- series::VecOrMat{T}: collection of vectors
- rep_value::Symbol: value to be used for representation of the aggregated interval (we use mean by default)
"""
# (second method)
function aggreg1D(series::VecOrMat{T}, representation::Symbol = :mean) where {T<:Float64}
    #### Function to aggregate vectors under a representative value, normally used in hierarchical
    #### clustering. It merges a set of vectors (series) merging the state (i) with (i + 1) whenever there's a marker

    @assert length(size(series)) <= 2 "Series to be aggregated need to be represented as a Matrix (2D)"

    # Initialization
    lseries = size(series,1)                      # Length of vectors
    nseries = size(series,2)                      # Number of vectors
    D = 1:nseries                                 # Range for the number of series
    class = collect(1:lseries)                    # Vector to mark the new classes
    new_v = copy(series)                          # New vector to be formed after aggregation

    # Calculate the representative value and replace it in the new_v
    if representation == :mean
        for d in D
            new_v[:,d] .= mean(series[:,d])
        end
    elseif representation == :medoid
        for d in D
            dist = pairwise(Euclidean(),series[:,d])
            s_medoid = kmedoids(dist, 1)
            s_medoid = s_medoid.medoids[1]
            new_v[:,d] .= series[:,d][s_medoid]
        end
    else
        @assert false "Representation method not defined."
    end

    return new_v, class
end

"""
aggreg1D(series::VecOrMat{T}, marker::Vector{Bool}, rep_value::Symbol = :mean)
The function aggregates N one-dimensional vectors concomitantly following a marker gave as the input of the function. It assumes the collection of 1D vectors is vertically concatanated so forming a 2D matrix with the first dimension equal to the size of the vectors analysed and the second dimension equal to the number of vectors analysed.
Fields:
- series::VecOrMat{T}: collection of vectors
- marker::Vector{Bool}: vector indicating the assignment of the series
- rep_value::Symbol: value to be used for representation of the aggregated interval (we use mean by default)
"""
# (third method)
function aggreg1D(series::VecOrMat{T}, marker::Vector{Int}, representation::Symbol = :mean) where {T<:Float64}
    #### Function to aggregate vectors under a representative value, normally used in hierarchical
    #### clustering. It merges a set of vectors (series) merging the state (i) with (i + 1) whenever they are equal

    @assert length(size(series)) <= 2 "Series to be aggregated need to be represented as a Matrix (2D)"

    # Initialization
    lseries = size(series,1)                      # Length of vectors
    nseries = size(series,2)                      # Number of vectors
    D = 1:nseries                                 # Range for the number of series
    new_lseries = maximum(marker)                 # New vectors length

    # Error control
    @assert lseries == length(marker) "Different sizes between series and assignment vector"

    # Finding the unique class vector
    class = unique(sort(marker))
    weights = [count(i->(i == j),marker) for j in class] .|> Int

    # Declaring the new VecOrMat object
    new_v = zeros(length(class),nseries)                          # New vector to be formed after aggregation

    # Calculate the representative value and replace it in the new_v
    if representation == :mean
        for i in 1:length(class), d in D
            new_v[i,d] = mean(series[marker .== class[i],d])
        end
    elseif representation == :medoid
        for i in 1:class[end], d in D
            dist = pairwise(Euclidean(),series[marker .== class[i],d])
            s_medoid = kmedoids(dist, 1)
            s_medoid = s_medoid.medoids[1]
            new_v[i,d] = series[marker .== class[i],d][s_medoid]
        end
    else
        @assert false "Representation method not defined."
    end

    return new_v, class, weights
end

"""
cdad(u::Vector{<:Float64},v::Vector{<:Float64})
Function to compute the cumulative distribution absolute difference from vector `u`
to vector `v` assuming both have an equal mass (i.e., both normalised) and
they are ordered.
"""
function cdad(u::Vector{T},v::Vector{T}) where {T<:Float64}
    # Error control
    @assert length(u) == length(v) "Vectors need to have an equal length."
    # @assert round(sum(u);digits=3) ≈ round(sum(v);digits=3) "Vectors need to have an equal mass."
    
    # Compute the differences between pairs of successive values from u and v.
    deltas = u .- v

    # Calculate the CDFs of u and v using their weights, if specified.
    cdf = sum(abs.(cumsum(deltas)))

    return cdf
end

"""
load_series_instance(series::VecOrMat{T}, block_size::S, current_k::S, stopping_k::S = 1, dm::Symbol = :ed, rep_value::Symbol = :mean, lseries::S = size(series,1), nseries::S = size(series,2)
) where {T<:Float64, S<:Int}
Function to load the data provided and fit into a struct SeriesInstance().
"""
function load_series_instance(series::VecOrMat{T},
    block_size::S,
    current_k::S,
    stopping_k::S = 1,
    dm::Symbol = :ed,
    rep_value::Symbol = :mean,
    lseries::S = size(series,1),
    nseries::S = size(series,2),
    series_dc::VecOrMat{T} = copy(series),
    ord_dc::VecOrMat{S} = repeat(collect(1:lseries),1,nseries)
) where {T<:Float64, S<:Int}

    @assert block_size <= current_k - stopping_k + 1 "Aggregation not possible: block_size $block_size and stopping_k $stopping_k need to be checked."

    # Create an object of type SeriesInstance
    _SeriesInstance = SeriesInstance(series,block_size,stopping_k,current_k,dm,rep_value,lseries,nseries,series_dc,ord_dc)

    return _SeriesInstance
end

"""
load_clust_instance(k_cent::VecOrMat{T}, weights::Vector{S} = ones(size(k_cent,1)), series_clust::Vector{S}, search_range::UnitRange = collect(1:size(k_cent,1))
) where {T<:Float64,S<:Int}
Function to load the clusters attributes.
"""
function load_clust_instance(k_cent::VecOrMat{T},
    series_clust::Vector{S},
    weights::Vector{S} = ones(size(k_cent,1)),
    search_range::UnitRange = collect(1:size(k_cent,1)),
    dc_mode::Bool = false
) where {T<:Float64,S<:Int}

    nclusters = size(k_cent,1)

    # Create an object of type ClustInstance
    _ClustInstance = ClustInstance(k_cent,weights,series_clust,nclusters,search_range, dc_mode)

    return _ClustInstance
end

"""
sorting_custom(series::VecOrMat{T},merging_interval::UnitRange,series_dc::VecOrMat{T},ord_dc::VecOrMat{S}
) where {T<:Number,S<:Int}

Customised function to replace `sort()` when the initial series_dc is known in advance.
"""
function sorting_custom(series::VecOrMat{T},merging_interval::UnitRange,series_dc::VecOrMat{T},ord_dc::VecOrMat{S}
) where {T<:Number,S<:Int}
    ## TODO: Implement methodology for mediods or other representative points

    # Declaring ordered series without merged values
    short_series_dc = Matrix{Float64}(undef,size(series,1)-length(merging_interval),size(series,2))

    # Return object
    ret_series_dc = copy(series_dc)

    for j in 1:size(series,2)
        # Forming plateaus
        k_mean = mean(series[merging_interval,j])
        
        # Forming new series_dc without merged values
        short_series_dc[:,j] = series[ord_dc[.![i in merging_interval for i in ord_dc[:,j]],j],j]

        # Finding where to allocate the new k_mean
        line_aux = 1 |> Int       # Number of the line tested
        keep_it = true            # Iteration auxiliar variable

        # Forming ordered series
        while keep_it && line_aux < size(short_series_dc,1)
            # If we reached the point where the centroid should be ordered or the end and we haven't replaced any lines we add a repetition of
            # k_mean equal to the size of the merging_interval
            if short_series_dc[line_aux,j] <= k_mean || line_aux == size(short_series_dc,1)
                ret_series_dc[line_aux:end,j] = vcat(repeat([k_mean],length(merging_interval)),short_series_dc[line_aux:end,j])
                keep_it = false
            else
                ret_series_dc[line_aux,j] = short_series_dc[line_aux,j]
            end

            # Keep iterating
            line_aux += 1
        end    
    end

    return ret_series_dc
end

"""
search_min_dist(_SeriesInstance, _ClustInstance)
Search minimal distance and return the marker, the min_dist position, and the clusters to be merged.
"""
function search_min_dist(_SeriesInstance, _ClustInstance)
    # Unpacking SeriesInstance
    series = _SeriesInstance.series
    block_size = _SeriesInstance.block_size
    dm = _SeriesInstance.dm
    rep_value = _SeriesInstance.rep_value
    lseries = _SeriesInstance.lseries
    nseries = _SeriesInstance.nseries
    series_dc = _SeriesInstance.series_dc
    ord_dc = _SeriesInstance.ord_dc

    # Unpacking ClustInstance
    series_clust = _ClustInstance.series_clust
    search_range = _ClustInstance.search_range
    dc_mode = _ClustInstance.dc_mode

    # Vector with distances (to be updated as it goes)
    dist = Vector{Float64}(undef,length(search_range))

    # Sets
    N = 1:nseries
    K = copy(search_range)

    # Compute the distance for each aggregation (i.e., changing the merging_clust)
    # TODO: implement parallelisation such as '@async Threads.@threads @inbounds for k in K'
    @inbounds Threads.@threads for k_search in K
        # Merging to be tested (neighbouring hypothesis)
        # TODO: implement non-neighbouring hypothesis
        # merging_clust = k:k+block_size-1
        # Create a temporary marker to merge the clusters tested
        marker_temp = [sc in k_search:k_search+block_size-1 for sc in series_clust]

        # Separation needed for duration curves analysis
        if dc_mode
            # Using duration curves as the comparison series
            series_comp = copy(series_dc)
            # Forming the centroids in a decrescent order
            k_cent_comp = sorting_custom(series,minimum(collect(1:lseries)[marker_temp]):maximum(collect(1:lseries)[marker_temp]),series_dc,ord_dc)
        else
            # Part of series compared
            series_comp = series[marker_temp,:]
            # Centroids of the temporarily formed cluster (TODO: implement another method for aggreg1D receiving the clusters with respective weights)
            (k_cent_comp,) = aggreg1D(series_comp, rep_value)
        end

        # Distance computation
        dist[k_search] = compute_dist(N, dm, series_comp, k_cent_comp)
    end

    # Find whenever the min_dist occurs first (i.e., using findmin()[2])
    ## TODO: implement the multiple merges (e.g., using findall())
    min_dist = findmin(dist)[2] |> Int
    merging_clust = min_dist:min_dist+block_size-1

    # Create a flag to the positions in the series that will be aggregated
    marker = update_marker(_SeriesInstance, _ClustInstance, min_dist)

    return marker, min_dist, merging_clust
end

"""
find_clusters!(_SeriesInstance, _ClustInstance, _DistUpdate)
Function calling the necessary functions for the data aggregation procedure (to be enveloped in a loop with the necessary conditions, e.g., while nclusters > stopping_k).
"""
function find_clusters!(_SeriesInstance, _ClustInstance, _DistUpdate)

    # Find the min distance and create a marker
    (marker, min_dist, merging_clust) = search_min_dist(_SeriesInstance, _ClustInstance)

    # Update _DistUpdate dictionary with the minimal distance found and the new marker
    _DistUpdate = merge(+, _DistUpdate, Dict(marker => DistUpdate(min_dist,merging_clust)))

    # Update clusters and series_clust
    update_clust!(_ClustInstance, _SeriesInstance, min_dist)

    # Update number of clusters k
    new_current_k = _ClustInstance.nclusters
    update_k!(_SeriesInstance, new_current_k)

    return new_current_k, _DistUpdate
end

"""
function compute_dist(N::UnitRange, dm::Symbol, series_comp::VecOrMat{T}, k_cent_comp::VecOrMat{T})  where {T <: Float64}
Compute the distances between two matrices (series_comp and k_cent_comp).
"""
function compute_dist(N::UnitRange, dm::Symbol, series_comp::VecOrMat{T}, k_cent_comp::VecOrMat{T})  where {T <: Float64}
    # Distances matrix designation (TODO: add more discrepancy metrics (e.g., DTW and others) and wrap their calculation in a function)
    if dm == :wd  # Wasserstein distance (or cumulative distribution absolute difference) between the agglomerated cluster and the original series
        dist = sum(cdad(series_comp[:,n], k_cent_comp[:,n]) for n in N)
    elseif dm == :ed  # Ward's cluster criterion (min variance proxy) for the original series when clustering `series_range`
        dist = sum(sum((series_comp[:,n] .- k_cent_comp[:,n]).^2 for n in N))
    else
        @assert false "Distance metric not defined."
    end

    return dist
end

"""
update_marker(_SeriesInstance, _ClustInstance, min_dist::T) where {T <: Int}
Update the clusters assignment (Bool marker indicating which positions are to be merged)
"""
## TODO: implement the case of simultaneous multiple merges
function update_marker(_SeriesInstance, _ClustInstance, min_dist::T) where {T <: Int}
    block_size = _SeriesInstance.block_size
    lseries = _SeriesInstance.lseries
    series_clust = _ClustInstance.series_clust

    # Create the new marker with 1 for the elements within the series to be aggregated
    marker = [series_clust[i] in min_dist:min_dist+block_size-1 for i in 1:lseries] |> Vector{Bool}
    
    return marker
end

## TODO: normalisation function (for example to use medoids representative values and the wasserstein distance as the discrepancy metric - via cdad)

"""
function replace_lines(mat::VecOrMat{T}, mergings::OrderedDict{UnitRange{Int64},Matrix{T}}) where {T <: Float64}
Replace lines of a matrix/vector by another matrix/vector with the same number of columns than the original matrix/vector.
"""
function replace_lines(mat::VecOrMat{T}, mergings::OrderedDict{UnitRange{Int64},Matrix{T}}) where {T <: Float64}
    nmergings = length(mergings)
    nrows = size(mat,1)
    ncols = size(mat,2)
    merg_keys = keys(mergings) .|> UnitRange
    
    # Error control for overlapping intervals and matrix scope
    lmerged = [] |> Vector{Int}
    for i in 1:nmergings
        lmerged = vcat(lmerged,collect(merg_keys[i]))
    end
    lmerged = sort(lmerged)
    @assert maximum(count(i -> i == lmerged[j], lmerged) for j in 1:length(lmerged)) == 1 "Overlapping merging intervals"
    @assert all([i in 1:nrows for i in lmerged]) "Lines to be merged out of matrix scope"

    # Flag to indicate if the merging was already transfered for the final matrix
    flag_merg = zeros(nmergings) |> Vector{Bool}
    
    # New matrix
    new_mat = Array{Float64}(undef,0,ncols)

    # Loop to merge the new lines
    for i in 1:nrows
        merg_pos = findmax([i in merg_keys[j] for j in 1:nmergings])[2]       # find the position in `mergings` where this merging happens
        if (i in lmerged) && (flag_merg[merg_pos] == 0)                       # case in which the line is in the mergings list and still not added to the new matrix
            new_mat = vcat(new_mat,mergings[merg_keys[merg_pos]])
            flag_merg[merg_pos] = 1
        elseif !(i in lmerged)
            new_mat = vcat(new_mat,reshape(mat[i,:],(1,ncols)))
        end
    end
    
    return new_mat
end

"""
function update_clust!(_ClustInstance, _SeriesInstance, min_dist::T) where {T <: Int}
Updates clusters attributes in the _ClustInstance object.
"""
function update_clust!(_ClustInstance, _SeriesInstance, min_dist::T) where {T <: Int}
    series = _SeriesInstance.series
    series_clust = _ClustInstance.series_clust
    nclusters = _ClustInstance.nclusters
    k_cent = _ClustInstance.k_cent

    block_size = _SeriesInstance.block_size
    rep_value = _SeriesInstance.rep_value
    nseries = _SeriesInstance.nseries
    lseries = _SeriesInstance.lseries

    # Ranges
    L = 1:lseries
    merging_clust = min_dist:min_dist + block_size - 1

    @assert length(merging_clust) >= 2 "unitary merging_clust not allowed"

    # Updating the centroids
    (k_change,) = aggreg1D(series[[sc in merging_clust for sc in series_clust],:], rep_value)
    k_change = reshape(k_change[1,:],(1,nseries))
    k_cent = replace_lines(k_cent,sort(Dict(merging_clust => k_change)))

    # Update series_clust
    change_marker = zeros(lseries) |> Vector{Int}                        # Mark whenever the change in the cluster assignment happens for the series_clust
    @inbounds for l in 2:L[end]
        if (series_clust[l] in merging_clust[2:end]) && (series_clust[l] != series_clust[l-1])
            change_marker[l] = 1
        end
    end

    # Transform change_marker to carry the cumulative change in the clusters order
    c_change = 0                        # Change in the number of clusters
    @inbounds for l in L
        c_change = c_change - change_marker[l]
        change_marker[l] = c_change
    end
    
    series_clust = (series_clust .+ change_marker) |> Vector{Int64}
    series_clust = series_clust |> Vector{Int64}

    # New number of clusters
    nclusters = series_clust[end]
    # Update weights
    weights = [sum(series_clust .== i) for i in 1:nclusters] .|> Int
    # Update search_range
    new_search_range = 1:(1+nclusters-block_size)

    # Update struct
    _ClustInstance.series_clust = series_clust
    _ClustInstance.weights = weights
    _ClustInstance.nclusters = nclusters
    _ClustInstance.k_cent = k_cent
    _ClustInstance.search_range = new_search_range

    return
end

"""
update_k!(_SeriesInstance, new_current_k::Int)
Updates number of clusters in the _SeriesInstance object.
"""
function update_k!(_SeriesInstance, new_current_k::Int)
    _SeriesInstance.current_k = new_current_k
    stopping_k = _SeriesInstance.stopping_k

    @assert new_current_k >= stopping_k "Number of clusters $new_current_k is now < then $stopping_k"
end

"""
"""
function write_clust_instance!(steps_per_block::Int,num_hours::Int,rep::String,dm::String,instances_path::AbstractString, ParamsDict, ClustDict, SeriesDict, ClustersRange::Vector{Int}, num_nodes::Int)
    for i in ClustersRange
        num_clusters = i
        # Declaring the instance
        if dm == "ward"
            if rep == "mean"
                instance = string(lpad(num_nodes,2,"0"),"n",num_hours,"h",lpad(num_clusters,4,"0"),"c",steps_per_block,"b_","m","m")
            elseif rep == "medoid"
                instance = string(lpad(num_nodes,2,"0"),"n",num_hours,"h",lpad(num_clusters,4,"0"),"c",steps_per_block,"b_","m","d")
            end
        elseif dm == "wd"
            if rep == "mean"
                instance = string(lpad(num_nodes,2,"0"),"n",num_hours,"h",lpad(num_clusters,4,"0"),"c",steps_per_block,"b_","w","m")
            elseif rep == "medoid"
                instance = string(lpad(num_nodes,2,"0"),"n",num_hours,"h",lpad(num_clusters,4,"0"),"c",steps_per_block,"b_","w","d")
            end
        end
        clust_instance_path = joinpath(instances_path,instance)
        mkpath(clust_instance_path);

        # Creating the clusters instance
        _ClustInstance = ClustDict[string(num_clusters)]
        _SeriesInstance = SeriesDict[string(num_clusters)];

        # Writing the weights in a csv
        weights_df = DataFrame(Weights = _ClustInstance.weights) |> CSV.write(joinpath(clust_instance_path,"weights.csv"));

        # Defining the path to data related to n and t
        nodes_path = joinpath(instances_path,instance,"nodes")
        mkpath(nodes_path)
        cd(nodes_path)

        # FTR series
        series = _SeriesInstance.series

        # Writing csv's
        for n in 1:num_nodes
            Demand = _ClustInstance.k_cent[:,n]
            Avail_Wind_On = _ClustInstance.k_cent[:,num_nodes + n]
            Avail_Wind_Off = _ClustInstance.k_cent[:,2*num_nodes + n]
            Avail_Sol = _ClustInstance.k_cent[:,3*num_nodes + n]
            (Hyd_In,) = aggreg1D(ParamsDict["AH_nt"][n,:]|>Vector{Float64},_ClustInstance.series_clust)
            (HydRoR_In,) = aggreg1D(ParamsDict["AR_nt"][n,:]|>Vector{Float64},_ClustInstance.series_clust)
            node_df = DataFrame(Demand=Demand, Avail_Sol=Avail_Sol, Avail_Wind_On=Avail_Wind_On, Avail_Wind_Off=Avail_Wind_Off, Hyd_In=Hyd_In[:,1], HydRoR_In=HydRoR_In[:,1]) |> CSV.write("$n.csv")
        end
        
        # Writing representative periods .JSON
        rep_periods = Dict("T" => num_clusters)
        save_json(rep_periods,joinpath(instances_path,instance,"rep_periods.json"))
    end
end

"""
Read different clustering instances and return time-dependent parameters.
"""
function read_clust_instance(clust_method::AbstractString, clust_method_path::AbstractString, instance_clust::AbstractString, ParamsDict::Dict{String, Any}; nosun::Bool = false)
    # Unpack parameters needed from ParamsDict
    G = ParamsDict["G"]
    N = ParamsDict["N"]

    # Loading clusters dictionary
    ClustDict = load(clust_method_path)

    # Number of representative hours from the instance_clust
    rep_hours = parse(Int,instance_clust[9:12])

    if occursin("day",clust_method)
        # Number of days from the instance_clust
        hours_day = 24
        days_year = 365

        @assert rep_hours%hours_day == 0 "Non-integer number of days"
        num_days = rep_hours/hours_day |> Int
        @assert num_days < days_year "More days than possible in an year"

        # Declaring clustering instance dictionary
        _ClustInstance = ClustDict[string(num_days)]

        # Change the num_clusters from num_days to num_hours
        _ClustInstance.nclusters = hours_day*num_days
        # Declare rep. periods
        T = 1:_ClustInstance.nclusters

        # Availability generation series to be considered (1: Wind-on, 2: Wind-off, 3: Solar)
        A_range = 1:3

        # Transforming k_cent (rep. days) into D_nt and A_gnt (here additional columns, e.g., RMSE columns, will not be considered)
        D_nt = reshape(permutedims(_ClustInstance.k_cent[:,1:length(N)*24],[2,1]),(length(N),hours_day*num_days))
        A_gnt = reshape(permutedims(_ClustInstance.k_cent[:,length(N)*hours_day+1:length(N)*24+length(A_range)*length(N)*hours_day],[2,1]),
        (length(A_range),length(N),hours_day*num_days));
        A_gnt[A_gnt .< 0.001] .= 0

        # Increase size of A_gnt to include availability of non-renewable generation techs.
        A_gnt_new = zeros(length(G),size(A_gnt,2),size(A_gnt,3))
        A_gnt_new[A_range,:,:] = A_gnt
        A_gnt_new[length(A_range)+1:length(G),:,:] .= 1
        A_gnt = A_gnt_new

        # Updating clusters weights and series_clust
        τ_t = zeros(Int,hours_day*length(_ClustInstance.weights))

        for i in 1:num_days
            τ_t[hours_day*(i-1)+1:hours_day*i] .= _ClustInstance.weights[i]
        end

        # Demand and availability matrices (arranged in days and hours)
        AH_dn = permutedims(reshape(ParamsDict["AH_nt"],(Int(length(N)*hours_day),days_year)),[2,1]) |> Matrix{Float64};
        AR_dn = permutedims(reshape(ParamsDict["AR_nt"],(Int(length(N)*hours_day),days_year)),[2,1]) |> Matrix{Float64};

        # New AH_dn matrix to store aggregated values (in rep. days)
        AH_dn_new = zeros(num_days,length(N)*hours_day)
        AR_dn_new = zeros(num_days,length(N)*hours_day)

        # Aggregate values in rep. days
        for i in 1:length(N)*hours_day
            (AH_dn_new[:,i],) = aggreg1D(AH_dn[:,i]|>Vector{Float64},_ClustInstance.series_clust)
            (AR_dn_new[:,i],) = aggreg1D(AR_dn[:,i]|>Vector{Float64},_ClustInstance.series_clust)
        end

        # Turn rep. days in rep. hours
        AH_nt = reshape(permutedims(AH_dn_new[:,1:length(N)*24],[2,1]),(length(N),hours_day*num_days))
        AR_nt = reshape(permutedims(AR_dn_new[:,1:length(N)*24],[2,1]),(length(N),hours_day*num_days))
    else
        # Declaring clustering instance dictionary
        _ClustInstance = ClustDict[string(rep_hours)]

        # Clusters weight
        τ_t = _ClustInstance.weights

        # Declare rep. periods
        # rep_periods = Dict("T" => num_clusters)
        T = 1:_ClustInstance.nclusters

        # Clusters weight
        τ_t = _ClustInstance.weights

        # Optimisation objects
        D_nt = zeros(length(N), length(T))
        A_gnt = ones(length(G), length(N), length(T))
        AH_nt = zeros(length(N), length(T))
        AR_nt = zeros(length(N), length(T))

        for n in N
            # Time-dependent parameters
            (D_nt[n, :],) = aggreg1D(ParamsDict["D_nt"][n,:]|>Vector{Float64},_ClustInstance.series_clust)
            (A_gnt[1, n, :],) = aggreg1D(ParamsDict["A_gnt"][1,n,:]|>Vector{Float64},_ClustInstance.series_clust)
            (A_gnt[2, n, :],) = aggreg1D(ParamsDict["A_gnt"][2,n,:]|>Vector{Float64},_ClustInstance.series_clust)
            (A_gnt[3, n, :],) = aggreg1D(ParamsDict["A_gnt"][3,n,:]|>Vector{Float64},_ClustInstance.series_clust)
            (AH_nt[n,:],) = aggreg1D(ParamsDict["AH_nt"][n,:]|>Vector{Float64},_ClustInstance.series_clust)
            (AR_nt[n,:],) = aggreg1D(ParamsDict["AR_nt"][n,:]|>Vector{Float64},_ClustInstance.series_clust)
        end
        # Rounding to improve numerical stability
        D_nt = round.(D_nt; digits = 5)
        A_gnt = round.(A_gnt;digits = 5)
        A_gnt[A_gnt .< 0.001] .= 0
        AH_nt = round.(AH_nt;digits=0)
        AR_nt = round.(AR_nt;digits=0)
    end

    # Producing nosun
    if nosun
        gsun = 3
        A_gnt[gsun,:,:] .= 0
    end

    # Forming parameters struct with updated T, τ_t, A_gnt, and D_nt
    parameters = Params(
        ParamsDict["region_n"], ParamsDict["max_dem_n"], ParamsDict["technology_g"], ParamsDict["G"], ParamsDict["G_r"], ParamsDict["N"], ParamsDict["L"], ParamsDict["L_ind"], 
        T, ParamsDict["S"], ParamsDict["H"], ParamsDict["κ"], ParamsDict["μ"], ParamsDict["C"], ParamsDict["C̄"], ParamsDict["C_E"], ParamsDict["R_E"], τ_t,
        ParamsDict["Gmin_gn"], ParamsDict["Gmax_gn"], A_gnt, D_nt, ParamsDict["I_g"], ParamsDict["M_g"], ParamsDict["C_g"], ParamsDict["e_g"], ParamsDict["E_g"],
        ParamsDict["r⁻_g"], ParamsDict["r⁺_g"], ParamsDict["I_l"], ParamsDict["M_l"], ParamsDict["C_l"], ParamsDict["B_l"], ParamsDict["e_l"], ParamsDict["Tmin_l"], ParamsDict["Tmax_l"],
        ParamsDict["ξ_s"], ParamsDict["I_s"], ParamsDict["C_s"], ParamsDict["Smin_sn"], ParamsDict["Smax_sn"], ParamsDict["Wmax_hn"], ParamsDict["Wmin_hn"], ParamsDict["Hmax_hn"],
        ParamsDict["Hmin_hn"], ParamsDict["HRmax_n"], ParamsDict["Fmin_n"], AH_nt, AR_nt, ParamsDict["I_h"], ParamsDict["M_h"], ParamsDict["C_h"],
        ParamsDict["e_h"], ParamsDict["E_h"], ParamsDict["r⁻_h"], ParamsDict["r⁺_h"]
    );
    return parameters
end