using Statistics, OrderedCollections, Dates, Parameters

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
"""
mutable struct ClustInstance{T<:Float64,S<:Int}
    k_cent::VecOrMat{T}
    weights::Vector{S}
    series_clust::Vector{S}
    nclusters::S
    search_range::UnitRange
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
    @assert sum(u) ≈ sum(v) "Vectors need to have an equal mass."
    
    # Compute the differences between pairs of successive values from u and v.
    deltas = u .- v

    # Calculate the CDFs of u and v using their weights, if specified.
    cdf = sum(abs.(cumsum(deltas)))

    return cdf
end

"""
load_series_instance(series::VecOrMat{T}, block_size::S, current_k::S, stopping_k::S = 1, dm::Symbol = :ward, rep_value::Symbol = :mean, lseries::S = size(series,1), nseries::S = size(series,2)
) where {T<:Float64, S<:Int}
Function to load the data provided and fit into a struct SeriesInstance().
"""
function load_series_instance(series::VecOrMat{T},
    block_size::S,
    current_k::S,
    stopping_k::S = 1,
    dm::Symbol = :ward,
    rep_value::Symbol = :mean,
    lseries::S = size(series,1),
    nseries::S = size(series,2)
) where {T<:Float64, S<:Int}

    @assert block_size <= current_k - stopping_k + 1 "Aggregation not possible: block_size $block_size and stopping_k $stopping_k need to be checked."

    # Create an object of type SeriesInstance
    _SeriesInstance = SeriesInstance(series,block_size,stopping_k,current_k,dm,rep_value,lseries,nseries)

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
    search_range::UnitRange = collect(1:size(k_cent,1))
) where {T<:Float64,S<:Int}

    nclusters = size(k_cent,1)

    # Create an object of type SeriesInstance
    _ClustInstance = ClustInstance(k_cent,weights,series_clust,nclusters,search_range)

    return _ClustInstance
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
    nseries = _SeriesInstance.nseries

    # Unpacking ClustInstance
    series_clust = _ClustInstance.series_clust
    search_range = _ClustInstance.search_range

    # Vector with distances (to be updated as it goes)
    dist = Vector{Float64}(undef,length(search_range))

    # Sets
    N = 1:nseries
    K = copy(search_range)

    # Compute the distance for each aggregation (i.e., changing the merging_clust)
    @inbounds for k in K
        # Merging to be tested
        merging_clust = k:k+block_size-1
        # Create a temporary marker to merge the clusters tested
        marker_temp = [sc in merging_clust for sc in series_clust]
        # Part of series compared
        series_comp = series[marker_temp,:]
        # Centroids of the temporarily formed cluster (TODO: implement another method for aggreg1D receiving the clusters with respective weights)
        (k_cent_comp,) = aggreg1D(series_comp, rep_value)

        # Distance computation
        dist[k] = compute_dist(N, dm, series_comp, k_cent_comp)
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
    elseif dm == :ward  # Ward's cluster criterion (min variance proxy) for the original series when clustering `series_range`
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
update_k!(_SeriesInstance, new_current_k::Int)
Updates number of clusters in the _SeriesInstance object.
"""
# function generate_instance(parameters::Params, ClustDict, SeriesDict, clustering_levels::Vecto{Int}, output_path)