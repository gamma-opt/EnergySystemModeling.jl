using EnergySystemModeling
using FileIO
using JLD2
using Logging
using Dates

WRKDIR = "/scratch/work/condeil1/EnergySystemModeling.jl"

@info "Creating aggreg_TS directory"
output_dir = "aggreg_out"
mkpath(output_dir)

@info "Loading parameters"
examples_dir = joinpath(WRKDIR,"examples")
constants_path = joinpath(examples_dir,"constants")
structures_path = "8nodes"
instance = "big"
instances_path = joinpath(examples_dir,"structures",structures_path,"instances",instance)

parameters = Params(constants_path, instances_path)

@info "Declaring instances"
# Number of availability series to be considered in the clustering:
navail = 3

# Normalising demand series
dnt = transpose(parameters.D_nt)./repeat(transpose(maximum(parameters.D_nt[:,h] for h in 1:size(transpose(parameters.D_nt),1))),size(transpose(parameters.D_nt),1)) |> Array{Float64};

# Concatenating series
series = hcat(dnt, reshape(permutedims(parameters.A_gnt[1:navail,:,:],[3,1,2]),(size(dnt,1),navail*size(dnt,2)))) |> Array{Float64};

## Series attributes:
block_size = 2
stopping_k = 8750
current_k = size(series,1)
dm = :wd
rep_value = :mean
lseries = size(series,1)
nseries = size(series,2);

## Clustering attributes:
k_cent = copy(series)
weights = ones(lseries) |> Vector{Int}
series_clust = collect(1:lseries)
nclusters = lseries
search_range = 1:(1+lseries-block_size);

# Declare instances
_SeriesInstance = load_series_instance(series,
    block_size,
    current_k,
    stopping_k,
    dm,
    rep_value,
    lseries,
    nseries
)

_ClustInstance = load_clust_instance(k_cent,
    series_clust,
    weights,
    search_range
)

# Define a new copy method to copy _ClustInstance and _SeriesInstance
Base.copy(x::T) where T = T([getfield(x, k) for k ∈ fieldnames(T)]...)

# Initialize:
k = parameters.T[end]

# Dictionary to keep the min distances and respective markers/min_dist found in each iteration
_DistUpdate = Dict{Vector{Bool}, DistUpdate}()

# Dictionaries to store series_clust and k_cent
_ClustUpdate = Dict{String,ClustInstance}()
_SeriesUpdate = Dict{String,SeriesInstance}()

@info "Started clustering..."
while k >= stopping_k + block_size - 1
    # Intermediate saving and logging
    if k % 200 == 0
        @info string("Last update: ", round(100*k/parameters.T[end], digits=2), "% @ ", now())
        save(joinpath(output_dir,"clust_prelim.jld2"),_ClustUpdate)
        save(joinpath(output_dir,"series_prelim.jld2"),_SeriesUpdate)
    end
    global (k, _DistUpdate) = find_clusters!(_SeriesInstance, _ClustInstance, _DistUpdate)
    
    # Store series_clust and k_cent
    global _ClustUpdate = merge(+,_ClustUpdate,Dict("$k" => copy(_ClustInstance)))
    global _SeriesUpdate = merge(+,_SeriesUpdate,Dict("$k" => copy(_SeriesInstance)))
end

save(joinpath(output_dir,"clust_out.jld2"),_ClustUpdate)
save(joinpath(output_dir,"series_out.jld2"),_SeriesUpdate)