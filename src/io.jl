using CSV, JSON, DataFrames, JLD, MAT, FileIO, JLD2
include("aggreg.jl")

"""Equivalent annual cost (EAC)
# Arguments
* `cost::Real`: Net present value of a project.
* `n::Integer`: Number of payments.
* `r::Real`: Interest rate.
"""
function equivalent_annual_cost(cost::Real, n::Integer, r::Real)
    factor = if iszero(r) n else (1 - (1 + r)^(-n)) / r end
    return cost / factor
end

"""Loads parameter values for an instance from CSV and JSON files. Reads the following files from `instance_path`.
- `indices.json` with fields `G`, `G_r`, `N`, `L`, `T`, `S`
- `constants.json` with fields `kappa`, `C`, `C_bar`, `r`
- `nodes/` -- Time clustered data from the nodes with fields `Dem_Inc`, `Load_mod`, `Max_Load`, `Avail_Win`, `Avail_Sol`
  - `1.csv`
  - `2.csv`
  - ...
- `technology.csv` with fields `cost`, `lifetime`, `M`, `fuel_cost_1`, `fuel_cost_2`, `r_minus`, `r_plus`
- `transmission.csv` with fields `M`, `cost`, `dist`, `lifetime`, `C`, `B`
- `storage.csv` with fields `xi`, `cost`, `lifetime`, `C`, `b0_1, ..., b0_n`
# Arguments
- `DataInput_path::AbstractString`: Path to the instance directory.
"""
## Params for small instance

function Params(DataInput_path::AbstractString, Instances_path::AbstractString)
   
    # Load indexes and constant parameters
    indices = JSON.parsefile(joinpath(Instances_path, "indices.json"))

    # TODO: implement time period clustering: T, τ_t

    # Load indices. Convert JSON values to right types.
    G = indices["G"] |> Array{Int}
    G_r = indices["G_r"] |> Array{Int}
    N = indices["N"] |> Array{Int}
    L = indices["L"] |> Array{Array{Int}}
    L_ind = indices["L_ind"] |> Array{Int}
    T = 1:indices["T"]
    S = indices["S"] |> Array{Int}
    H = indices["H"] |> Array{Int}

    # Load constant parameters
    constants = JSON.parsefile(joinpath(DataInput_path, "constants.json"))
    κ = constants["kappa"]
    μ = constants["mu"]
    C = constants["C"]
    C̄ = constants["C_bar"]
    C_E = constants["C_E"]
    interest_rate = constants["r"]
    R_E = constants["R_E"]

    # Load time clustered parameters and region
    D_nt = zeros(length(N), length(T))
    A_gnt = ones(length(G), length(N), length(T))
    AH_nt = zeros(length(N), length(T))
    AR_nt = zeros(length(N), length(T))
    region_n = Array{AbstractString, 1}(undef, length(N))    

    # Load generation parameters (per node / per generation technology)
    Gmin_gn = zeros(length(G), length(N))
    Gmax_gn = zeros(length(G), length(N))
    gen_capacity = CSV.File(joinpath(Instances_path, "gen_capacity.csv")) |> DataFrame

    for n in N
        # Load node values from CSV files.
        nodes = CSV.File(joinpath(Instances_path, "nodes", "$n.csv")) |> DataFrame
        for g in G
            line_search = findall((gen_capacity.gen_tech .== g) .& (gen_capacity.node .== n))[1]
            Gmin_gn[g,n] = gen_capacity.gcap_min[line_search]
            Gmax_gn[g,n] = gen_capacity.gcap_max[line_search]
        end
        D_nt[n, :] = nodes.Demand[T]
        A_gnt[1, n, :] = nodes.Avail_Wind_On[T]
        A_gnt[2, n, :] = nodes.Avail_Wind_Off[T]
        A_gnt[3, n, :] = nodes.Avail_Sol[T]
        A_gnt[A_gnt .< 0.001] .= 0
        AH_nt[n,:] = nodes.Hyd_In[T]
        AR_nt[n,:] = nodes.HydRoR_In[T]
        region_n[n] =  nodes.Name[1]
    end

    # Load weights
    clust_weights = CSV.File(joinpath(Instances_path, "weights.csv")) |> DataFrame
    τ_t = clust_weights.Weights[T]

    # Load technology parameters
    gen_technology = joinpath(Instances_path, "gen_technology.csv") |>
        CSV.File |> DataFrame
    I_g = equivalent_annual_cost.(gen_technology.investment_cost .* 1000, gen_technology.lifetime,
                                  interest_rate) |> Array{AbstractFloat, 1}
    M_g = gen_technology.fixedOM .* 1000 |> Array{AbstractFloat, 1}
    C_g = gen_technology.fuel_cost ./ gen_technology.efficiency .+ gen_technology.varOM |> Array{AbstractFloat, 1}
    e_g = gen_technology.efficiency |> Array{AbstractFloat, 1}
    E_g = gen_technology.emissions |> Array{AbstractFloat, 1}
    r⁻_g = gen_technology.r_minus |> Array{AbstractFloat, 1}
    r⁺_g = gen_technology.r_plus |> Array{AbstractFloat, 1}
    technology_g = gen_technology.name

    # Load transmission parameters
    I_l = zeros(length(L_ind))
    M_l = zeros(length(L_ind))
    C_l = zeros(length(L_ind))
    B_l = zeros(length(L_ind))
    e_l = zeros(length(L_ind))
    Tmin_l = zeros(length(L_ind))
    Tmax_l = zeros(length(L_ind))
    transmission = joinpath(Instances_path, "transmission.csv") |>
        CSV.File |> DataFrame
    I_l = equivalent_annual_cost.(transmission.cost[1] .* transmission.dist .+ transmission.converter_cost[1],
                                  transmission.lifetime[1], interest_rate) |> Array{AbstractFloat, 1}
    M_l = transmission.M .* I_l |> Array{AbstractFloat, 1}
    C_l = transmission.C |> Array{AbstractFloat, 1}
    B_l = transmission.B |> Array{AbstractFloat, 1}
    e_l = transmission.efficiency |> Array{AbstractFloat, 1}
    Tmin_l = transmission.tcap_min |> Array{AbstractFloat, 1}
    Tmax_l = transmission.tcap_max |> Array{AbstractFloat, 1}

    # Load storage parameters
    ξ_s = zeros(length(S))
    I_s = zeros(length(S))
    C_s = zeros(length(S))
    Smin_sn = zeros(length(S),length(N))
    Smax_sn = zeros(length(S),length(N))
    storage = joinpath(Instances_path, "storage.csv") |>
        CSV.File |> DataFrame
    sto_capacity = joinpath(Instances_path, "sto_capacity.csv") |>
        CSV.File |> DataFrame
    for s in S
        ξ_s = storage.xi |> Array{AbstractFloat, 1}
        I_s = equivalent_annual_cost.(storage.cost .* 1000, storage.lifetime, interest_rate) |> Array{AbstractFloat, 1}
        C_s = storage.C |> Array{AbstractFloat, 1}
        for n in N
            line_search = findfirst((sto_capacity.s .== s) .& (sto_capacity.node .== n))
            Smin_sn[s,n] = sto_capacity.scap_min[line_search]
            Smax_sn[s,n] = sto_capacity.scap_max[line_search]
        end
    end
    Smin_sn = Smin_sn  |> Array{AbstractFloat, 2}
    Smax_sn = Smax_sn  |> Array{AbstractFloat, 2}

    # Load hydro capacity and technology parameters
    Wmax_n = zeros(length(N))
    Wmin_n = zeros(length(N))
    Hmin_n = zeros(length(N))
    Hmax_n = zeros(length(N))
    Fmin_n = zeros(length(N))
    HRcap_n = zeros(length(N))
    hydro = joinpath(Instances_path, "hydro.csv") |> CSV.File |> DataFrame;   
    Hmax_n[1:length(N)] = hydro.hcap_max[1:length(N)] |> Array{AbstractFloat, 1}
    Hmin_n[1:length(N)] = hydro.hcap_min[1:length(N)] |> Array{AbstractFloat, 1}
    HRcap_n[1:length(N)] = hydro.HydroRoR[1:length(N)] |> Array{AbstractFloat, 1}
    Wmax_n[1:length(N)] = hydro.wcap_max[1:length(N)] |> Array{AbstractFloat, 1}
    Wmin_n[1:length(N)] = hydro.wcap_min[1:length(N)] |> Array{AbstractFloat, 1}
    Fmin_n[1:length(N)] = hydro.Fmin[1:length(N)] |> Array{AbstractFloat, 1}
    hydro_technology = joinpath(Instances_path, "hydro_technology.csv") |> CSV.File |> DataFrame;   
    I_h = equivalent_annual_cost.(hydro_technology.investment_cost .* 1000, hydro_technology.lifetime,
                                    interest_rate) |> Array{AbstractFloat, 1}
    M_h = hydro_technology.fixedOM .* 1000 |> Array{AbstractFloat, 1}
    C_h = hydro_technology.fuel_cost ./ hydro_technology.efficiency .+ hydro_technology.varOM |> Array{AbstractFloat, 1}
    e_h = hydro_technology.efficiency |> Array{AbstractFloat, 1}
    E_h = hydro_technology.emissions |> Array{AbstractFloat, 1}
    r⁻_h = hydro_technology.r_minus |> Array{AbstractFloat, 1}
    r⁺_h = hydro_technology.r_plus |> Array{AbstractFloat, 1}
    
    # Return Params struct
    Params(
        region_n, technology_g, G, G_r, N, L, L_ind, T, S, κ, μ, C, C̄, C_E, R_E, τ_t, Gmin_gn, Gmax_gn, A_gnt, D_nt, I_g, M_g, C_g,
        e_g, E_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, e_l, Tmin_l, Tmax_l, ξ_s, I_s, C_s, Smin_sn, Smax_sn,
        Wmax_n, Wmin_n, Hmax_n, Hmin_n, HRcap_n, Fmin_n, AH_nt, AR_nt,
        I_h, M_h, C_h, e_h, E_h, r⁻_h, r⁺_h)
end


"""Save object into JSON file.
# Arguments
- `object`
- `output_path::AbstractString`: Full filepath, e.g., `path.json`.
"""
function save_json(object, filepath::AbstractString)
    open(filepath, "w") do io
        JSON.print(io, object)
    end
end

flatten(x::Array{<:Array, 1}) = Iterators.flatten(x)|> collect|> flatten
flatten(x::Array{<:Number, 1}) = x

shape(x::Array{<:Array, 1}) = vcat(shape(first(x)), [length(x)])
shape(x::Array{<:Number, 1}) = [length(x)]

mdim(x::Array{<:Array, 1}) = reshape(flatten(x), shape(x)...)
mdim(x::Array{<:Number, 1}) = x

transform(x::Array, t::Type{Array{T, N}}) where T <: Number where N = mdim(x)
transform(x::Array, t::Type{Array{T, 1}}) where T <: Array = x
transform(x::Number, t::Type{T}) where T <: Number = x
transform(x::Array, t::Type{Array{AbstractString, 1}}) = x

function convert_type(::Type{Array{T, N}}) where T <: Number where N
    t = T
    for i in 1:N
        t = Array{t, 1}
    end
    return t
end
convert_type(::Type{Array{T, 1}}) where T <: Array = Array{T, 1}
convert_type(t::Type{T}) where T <: Number = t
convert_type(::Type{Array{AbstractString, 1}}) = Array{AbstractString, 1}

"""Load values to type from JSON file.
# Arguments
- `type`
- `filepath::AbstractString`
"""
function load_json(type, filepath::AbstractString)
    objects = JSON.parsefile(filepath)
    fields = []
    for (s, t) in zip(fieldnames(type), fieldtypes(type))
        push!(fields,
              objects[string(s)] |>
              v -> convert(convert_type(t), v) |>
              v -> transform(v, t))
    end
    type(fields...)
end

#Replace NaN's with 0
function replace_nans!(array::Array{Float64, N}) where N
    for i in eachindex(array)
        if isnan(array[i])
            array[i] = zero(1)
        end
    end
end

"""Reads wind, solar and hydro data produced by the GlobalEnergyGIS package into CSV node files.
# Arguments
- `DataInput_path::AbstractString`: Path to the instance directory.
- `era_year::AbstractString`: Year of the ERA5 data used in the GlobalEnergyGIS package.
- `era_year::AbstractString`: Name of the GIS region used.
"""
function create_nodedata(DataInput_path::AbstractString, era_year::AbstractString, gisregion::AbstractString)

    #Read files
    solarvars = matread(joinpath(DataInput_path, "GISdata_solar$(era_year)_$gisregion.mat"))
    windvars = matread(joinpath(DataInput_path, "GISdata_wind$(era_year)_$gisregion.mat"))
    hydrovars = matread(joinpath(DataInput_path, "GISdata_hydro_$gisregion.mat"))
    demandvars = load(joinpath(instance_path, "SyntheticDemand_$(gisregion)_$(era_year).jld"), "demand")

    #Solar

    capacity_cspplantA = solarvars["capacity_cspplantA"]
    capacity_cspplantB = solarvars["capacity_cspplantB"]
    capacity_pvplantA = solarvars["capacity_pvplantA"]
    capacity_pvplantB = solarvars["capacity_pvplantB"]
    capacity_pvrooftop = solarvars["capacity_pvrooftop"]

    #Number of nodes
    n = size(capacity_cspplantA, 1)

    CFtime_cspplantA = solarvars["CFtime_cspplantA"]
    CFtime_cspplantB = solarvars["CFtime_cspplantB"]
    CFtime_pvplantA = solarvars["CFtime_pvplantA"]
    CFtime_pvplantB = solarvars["CFtime_pvplantB"]
    CFtime_pvrooftop = solarvars["CFtime_pvrooftop"]

    #Replace NaNs from hourly data with 0's
    replace_nans!(CFtime_cspplantA)
    replace_nans!(CFtime_cspplantB)
    replace_nans!(CFtime_pvplantA)
    replace_nans!(CFtime_pvplantB)
    replace_nans!(CFtime_pvrooftop)

    #Preallocate arrays for results
    avail_sol_cspA = zeros(Float64, T, n)
    avail_sol_cspB = zeros(Float64, T, n)
    avail_sol_pvA = zeros(Float64, T, n)
    avail_sol_pvB = zeros(Float64, T, n)
    avail_sol_pvrooftop = zeros(Float64, T, n)

    #Calculate absolute availability values and sum up the different classes for each solar power type
    for i in 1:T
        avail_sol_cspA[i,:,:] = permutedims(sum(capacity_cspplantA .* CFtime_cspplantA[i,:,:], dims=2))
        avail_sol_cspB[i,:,:] = permutedims(sum(capacity_cspplantB .* CFtime_cspplantB[i,:,:], dims=2))
        avail_sol_pvA[i,:,:] = permutedims(sum(capacity_pvplantA .* CFtime_pvplantA[i,:,:], dims=2))
        avail_sol_pvB[i,:,:] = permutedims(sum(capacity_pvplantB .* CFtime_pvplantB[i,:,:], dims=2))
        avail_sol_pvrooftop[i,:,:] = permutedims(sum(capacity_pvrooftop .* CFtime_pvrooftop[i,:,:], dims=2))
    end
    #Calculate total solar capacity
    capacity_solar = permutedims(sum(capacity_cspplantA, dims=2)) + permutedims(sum(capacity_cspplantB, dims=2)) + permutedims(sum(capacity_pvplantA, dims=2)) +
                     permutedims(sum(capacity_pvplantB, dims=2))  + permutedims(sum(capacity_pvrooftop, dims=2))
    #Sum the different solar types and divide with total capacity to get relative availability values
    avail_sol = (avail_sol_cspA + avail_sol_cspB + avail_sol_pvA + avail_sol_pvB + avail_sol_pvrooftop) ./ capacity_solar

    #Wind

    capacity_offshore = windvars["capacity_offshore"]
    capacity_onshoreA = windvars["capacity_onshoreA"]
    capacity_onshoreB = windvars["capacity_onshoreB"]

    CFtime_windoffshore = windvars["CFtime_windoffshore"]
    CFtime_windonshoreA = windvars["CFtime_windonshoreA"]
    CFtime_windonshoreB = windvars["CFtime_windonshoreB"]

    #Replace NaNs from hourly data with 0's
    replace_nans!(CFtime_windoffshore)
    replace_nans!(CFtime_windonshoreA)
    replace_nans!(CFtime_windonshoreB)

    #Preallocate arrays for results
    avail_wind_offshore = zeros(Float64, T, n)
    avail_wind_onshoreA = zeros(Float64, T, n)
    avail_wind_onshoreB = zeros(Float64, T, n)

    #Calculate absolute availability values and sum up the different classes for each wind power type
    for i in 1:T
        avail_wind_offshore[i,:,:] = permutedims(sum(capacity_offshore .* CFtime_windoffshore[i,:,:], dims=2))
        avail_wind_onshoreA[i,:,:] = permutedims(sum(capacity_onshoreA .* CFtime_windonshoreA[i,:,:], dims=2))
        avail_wind_onshoreB[i,:,:] = permutedims(sum(capacity_onshoreB .* CFtime_windonshoreB[i,:,:], dims=2))
    end
    #Calculate total wind capacity
    capacity_wind_on = permutedims(sum(capacity_onshoreA, dims=2)) + permutedims(sum(capacity_onshoreB, dims=2))
    capacity_wind_off = permutedims(sum(capacity_offshore, dims=2))
    #Sum the different wind types and divide with total capacity to get relative availability values
    avail_wind_on = (avail_wind_onshoreA + avail_wind_onshoreB) ./ capacity_wind_on
    avail_wind_off = avail_wind_offshore ./ capacity_wind_off

    #Hydro
    
    existingcapac = permutedims(hydrovars["existingcapac"]) .* 1000
    #Monthly inflow
    existinginflow = permutedims(hydrovars["existinginflowcf"])
    #Replace NaNs from monthly data with 0's
    replace_nans!(existinginflow)
    #Preallocate array for intermediate result
    avail_inflow = zeros(Float64, T, n)

    #Turn monthly inflow data into hourly data
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    lasthours = 24 * cumsum(days)
    firsthours = [1; 1 .+ lasthours[1:end-1]]
    for m in 1:12
        for i = firsthours[m]:lasthours[m]
            avail_inflow[i,:] = existinginflow[m,:]
        end
    end

    #Percentage of inflow that is to reservoirs. TODO: Read from file
    reservoirp = [0 0.953 1 0 0.433 0.261 0 0.456 0.069 0.838 0.563]

    #Hydro capacities
    capacity_hyd = existingcapac .* reservoirp
    capacity_hydRoR = existingcapac - capacity_hyd
    hydrocapacity = [capacity_hyd; capacity_hydRoR]
    
    #Reservoir and RoR inflow
    hyd_in = avail_inflow .* capacity_hyd
    hydRoR_in = avail_inflow .* capacity_hydRoR

    #Write a CSV file for each node
    for i in 1:n
        nodedata = zeros(T, 6)
        nodedata[:,1] = demandvars[:,i]
        nodedata[:,2] = avail_sol[:,i]
        nodedata[:,3] = avail_wind_on[:,i]
        nodedata[:,4] = avail_wind_off[:,i]
        nodedata[:,5] = hyd_in[:,i]
        nodedata[:,6] = hydRoR_in[:,i]
        nodedata = convert(DataFrame, nodedata)
        rename!(nodedata, ["Demand", "Avail_Sol", "Avail_Wind_On", "Avail_Wind_Off", "Hyd_In", "HydRoR_In"])
        CSV.write(joinpath(DataInput_path, "nodes", "$i.csv"), nodedata)
    end
end

"""
read_clusters(ClustersDataPath::AbstractString, k::Int)
Read clusters from a julia dictionary (.jld2) with the structure ClustInstance.
# Output:
- _ClustInstance: Instance of clustering with weights, centroids, etc.
- _SeriesInstance: Instance of series used to form the clusters.
"""
function read_clusters(ClustersDataPath::AbstractString, k::Int)
    _ClustUpdate = load(joinpath(ClustersDataPath,"clust_out.jld2"))
    _SeriesUpdate = load(joinpath(ClustersDataPath,"series_out.jld2"));

    return _ClustUpdate[string(k)], _SeriesUpdate[string(k)] 
end

"""
unpack_clusters_features(_ClustInstance::ClustInstance)
# Output:
- k_cent: centroids.
- weights: respective centroids' weights.
- series_clust: mapping between series and clusters.
"""
function unpack_clusters_features(_ClustInstance::ClustInstance)
    @assert _ClustInstance.nclusters == _ClustInstance.series_clust[end] "Number of clusters inconsistent in the cluster instance"

    k_cent = _ClustInstance.k_cent
    weights = _ClustInstance.weights
    series_clust = _ClustInstance.series_clust

    return k_cent, weights, series_clust
end