# Aggreg.jl

## Inputdata

Contains time series data (i.e., respecting a temporal order) with the following attributes:
* region: optional information to specify the region data belongs to
* period_from: date(year,month,day,hour) from which the data is from
* period_to: date(year,month,day,hour) until which the data is
* resolution: resolution of the original series (e.g., hours, days, weeks, etc.)
* nseries: number of time series inputed (treated separately as if independent from each other)
* data: Dictionary with an entry for each attribute `[file name (attribute: e.g technology)]-[column name (node: e.g. location)]`, Each entry of the dictionary is a 2-dimensional `time-steps T x periods K`-Array holding the data
  

## SeriesInstance

Contains the series and respective attributes needed for performing the aggregation with the following attributes: 
* series: original time series
* block_size: number of clusters to be aggregated at once each iteration
* stopping_k: reduced number of clusters pre-established
* current_k: current number of clusters
* dm: discrepancy metric to be used to measure the difference between clusters to be aggregated
* rep_value: function to assign a representative value for each cluster (by default, the mean of values inside the cluster)
* lseries: length of series
* nseries: number of series

## ClustInstance

Constains clusters information.

* k_cent: matrix with the representative points for each cluster
* weigths: 1-dimensional `periods K`-Array with the absolute weight for each period. E.g. for a yearly resolution time series of 365 days, sum(weights) = 365
* series_clust: clusters assignment for the original series (vector with the from-to relationship between the original series and the clusters)
* nclusters: current number of clusters
* search_range: range in which the clusters can be merged, considering the amount of clusters and the block_size (generally equals to 1: (nclusters-block_size+1))

## AggregInstance:

Contains the aggregation instance temporary values.

* merging_ clust
* series_comp
* k_cent_comp

## aggreg1D

The function aggregates N one-dimensional vectors concomitantly following a marker gave as the input of the function. It assumes the collection of 1D vectors is vertically concatanated so forming a 2D matrix with the first dimension equal to the size of the vectors analysed and the second dimension equal to the number of vectors analysed.

The function can use the following input:
* series: A collection of vectors‚
* marker: A vector indicating the assignment of the series
* rep_value: Value to be used for representation of the aggregated interval (we use mean by default)


## cdad

Function to compute the cumulative distribution absolute difference from vector `u`
to vector `v` assuming both have an equal mass (i.e., both normalised) and
they are ordered.

## load_series_instance

Function to load the data provided and fit into a struct SeriesInstance().

## load_clust_instance

Function to load the clusters attributes.

## search_min_dist

Search minimal distance and return the marker, the min_dist position, and the clusters to be merged.

## find_clusters!

Function calling the necessary functions for the data aggregation procedure (to be enveloped in a loop with the necessary conditions, e.g., while nclusters > stopping_k).

## compute_dist!

Compute the distances between two matrices (series_comp and k_cent_comp).


## update_marker

Update the clusters assignment (Bool marker indicating which positions are to be merged).

## replace_lines

Replace lines of a matrix/vector by another matrix/vector with the same number of columns than the original matrix/vector.

## update_clust!

Updates clusters attributes in the _ClustInstance object.

## update_k! 

Updates number of clusters in the _SeriesInstance object.

