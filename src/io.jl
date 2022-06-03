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

"""Loads parameter values for an instance from CSV and JSON files. Reads the following files from `instances_path`.
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
function Params(DataInput_path::AbstractString, instances_path_ftr::AbstractString)
    # Load indexes and constant parameters
    indices = JSON.parsefile(joinpath(instances_path_ftr, "indices.json"))

    # Load indices. Convert JSON values to right types.
    G = indices["G"] |> Array{Int}
    G_r = indices["G_r"] |> Array{Int}
    N = indices["N"] |> Array{Int}
    L = indices["L"] |> Array{Array{Int}}
    L_ind = indices["L_ind"] |> Array{Int}
    S = indices["S"] |> Array{Int}
    H = indices["H"] |> Array{Int}
    T = 1:indices["T"] |> Array{Int}

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

    # Load generation parameters (per node / per generation technology)
    Gmin_gn = zeros(length(G), length(N))
    Gmax_gn = zeros(length(G), length(N))
    gen_capacity = CSV.File(joinpath(instances_path_ftr, "gen_capacity.csv")) |> DataFrame

    for n in N
        # Load node values from CSV files.
        nodes = CSV.File(joinpath(instances_path_ftr, "nodes", "$n.csv")) |> DataFrame
        for g in G
            line_search = findfirst((gen_capacity.gen_tech .== g) .& (gen_capacity.node .== n))
            Gmin_gn[g,n] = gen_capacity.gcap_min[line_search]
            Gmax_gn[g,n] = gen_capacity.gcap_max[line_search]
        end
        D_nt[n, :] = round.(nodes.Demand[T]; digits = 5)
        A_gnt[1, n, :] = round.(nodes.Avail_Wind_On[T]; digits = 5)
        A_gnt[2, n, :] = round.(nodes.Avail_Wind_Off[T]; digits = 5)
        A_gnt[3, n, :] = round.(nodes.Avail_Sol[T]; digits = 5)
        A_gnt[A_gnt .< 0.001] .= 0
        AH_nt[n,:] = round.(nodes.Hyd_In[T]; digits = 0)
        AR_nt[n,:] = round.(nodes.HydRoR_In[T]; digits = 0)
    end

    # Nodes specifications
    region_n = Array{AbstractString, 1}(undef, length(N))    
    max_dem_n = Array{Float64, 1}(undef, length(N))    

    nodes_specs = CSV.File(joinpath(instances_path_ftr, "nodes_specs.csv")) |> DataFrame

    for n in N
        region_n[n] =  nodes_specs.Name[n]
        max_dem_n[n] = nodes_specs.Max_Demand[n]
    end

    # Load weights
    clust_weights = CSV.File(joinpath(instances_path_ftr, "weights.csv")) |> DataFrame
    τ_t = clust_weights.Weights[T]

    # Load technology parameters
    gen_technology = joinpath(instances_path_ftr, "gen_technology.csv") |>
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
    transmission = joinpath(instances_path_ftr, "transmission.csv") |>
        CSV.File |> DataFrame
    I_l = equivalent_annual_cost.(transmission.cost[1] .* transmission.dist .+ transmission.converter_cost[1],
                                  transmission.lifetime[1]|>Int64, interest_rate) |> Array{AbstractFloat, 1}
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
    storage = joinpath(instances_path_ftr, "storage.csv") |>
        CSV.File |> DataFrame
    sto_capacity = joinpath(instances_path_ftr, "sto_capacity.csv") |>
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
    Wmax_hn = zeros(length(H),length(N))
    Wmin_hn = zeros(length(H),length(N))
    Hmin_hn = zeros(length(H),length(N))
    Hmax_hn = zeros(length(H),length(N))
    Fmin_n = zeros(length(N))
    HRmax_n = zeros(length(N))

    hydro_capacity = joinpath(instances_path_ftr, "hydro_capacity.csv") |> CSV.File |> DataFrame   
    for h in H, n in N
        line_search = findfirst((hydro_capacity.hydro_tech .== h) .& (hydro_capacity.node .== n))
        Hmin_hn[h,n] = hydro_capacity.hcap_min[line_search]
        Hmax_hn[h,n] = hydro_capacity.hcap_max[line_search]
        Wmin_hn[h,n] = hydro_capacity.wcap_min[line_search]
        Wmax_hn[h,n] = hydro_capacity.wcap_max[line_search]
    end

    hydro = joinpath(instances_path_ftr, "hydro.csv") |> CSV.File |> DataFrame   
    HRmax_n[1:length(N)] = hydro.HydroRoR[1:length(N)] |> Array{AbstractFloat, 1}
    Fmin_n[1:length(N)] = hydro.hyd_flow_min[1:length(N)] |> Array{AbstractFloat, 1}

    hydro_technology = joinpath(instances_path_ftr, "hydro_technology.csv") |> CSV.File |> DataFrame;   
    I_h = equivalent_annual_cost.(hydro_technology.investment_cost .* 1000, hydro_technology.lifetime,
                                    interest_rate) |> Vector{Float64}
    M_h = hydro_technology.fixedOM .* 1000 |> Vector{Float64}
    C_h = hydro_technology.fuel_cost ./ hydro_technology.efficiency .+ hydro_technology.varOM |> Vector{Float64}
    e_h = hydro_technology.efficiency |> Vector{Float64}
    E_h = hydro_technology.emissions |> Vector{Float64}
    r⁻_h = hydro_technology.r_minus |> Vector{Float64}
    r⁺_h = hydro_technology.r_plus |> Vector{Float64}
    
    # Return Params struct
    Params(
        region_n, max_dem_n, technology_g, G, G_r, N, L, L_ind, T, S, H, κ, μ, C, C̄, C_E, R_E, τ_t, Gmin_gn, Gmax_gn, A_gnt, D_nt, I_g, M_g, C_g,
        e_g, E_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, e_l, Tmin_l, Tmax_l, ξ_s, I_s, C_s, Smin_sn, Smax_sn,
        Wmax_hn, Wmin_hn, Hmax_hn, Hmin_hn, HRmax_n, Fmin_n, AH_nt, AR_nt,
        I_h, M_h, C_h, e_h, E_h, r⁻_h, r⁺_h)
end

"""Perform function Params() and save objects (i.e., A_gnt, AH_nt, etc.) as a dictionary.
# Arguments
- `constants_path::AbstractString`: Path to the constants JSON file
- `instances_path_ftr::AbstractString`: Path to the instance directory
- `output_path::AbstractString`: Path where the dictionary with the parameters is saved.
"""
function perform_Params(constants_path::AbstractString, instances_path_ftr::AbstractString, params_path::AbstractString)
    # Create parameters object
    parameters = Params(constants_path, instances_path_ftr)

    # Create parameters dictionary folder
    mkpath(params_path)

    # Saving parameters into a dictionary
    parameters_dict = Dict(fieldnames(typeof(parameters)) .|> String .=> [getfield(parameters,i) for i in fieldnames(typeof(parameters))])
    JLD2.save(joinpath(params_path, "parameters.jld2"), parameters_dict)
end 

"""Read time-dependent parameters (i.e., T, D_nt, A_gnt, AH_nt, AR_nt) and create a new parameters struct.
# Arguments
- `params_path_ftr::AbstractString`: Path to the FTR parameters
- `instance_path_clust::AbstractString`: Path to the instance folder.
"""
function change_time_parameters(params_path_ftr::AbstractString, instance_path_clust::AbstractString; nosun::Bool = false)
    # Reading FTR parameters (assuming perform_Params was ran and the dictionary parameters.jld2 is available)
    ParamsDict = load(joinpath(params_path_ftr,"parameters.jld2"))

    # Declare representative periods
    rep_periods = JSON.parsefile(joinpath(instance_path_clust,"rep_periods.json"))
    T = 1:rep_periods["T"]
    ParamsDict["T"] = T

    # Load weights
    clust_weights = CSV.File(joinpath(instance_path_clust, "weights.csv")) |> DataFrame
    τ_t = clust_weights.Weights[T]
    ParamsDict["τ_t"] = τ_t

    # Read time-dependent parameters
    D_nt = zeros(length(ParamsDict["N"]), length(T))
    A_gnt = ones(length(ParamsDict["G"]), length(ParamsDict["N"]), length(T))
    AH_nt = zeros(length(ParamsDict["N"]), length(T))
    AR_nt = zeros(length(ParamsDict["N"]), length(T))

    for n in ParamsDict["N"]
            # Load node values from CSV files.
            nodes = CSV.File(joinpath(instance_path_clust,"nodes","$n.csv")) |> DataFrame
            D_nt[n, :] = round.(nodes.Demand[T]; digits = 5)
            A_gnt[1, n, :] = round.(nodes.Avail_Wind_On[T]; digits = 5)
            A_gnt[2, n, :] = round.(nodes.Avail_Wind_Off[T]; digits = 5)
            A_gnt[3, n, :] = round.(nodes.Avail_Sol[T]; digits = 5)
            A_gnt[A_gnt .< 0.001] .= 0
            AH_nt[n,:] = round.(nodes.Hyd_In[T]; digits = 0)
            AR_nt[n,:] = round.(nodes.HydRoR_In[T]; digits = 0)
    end

    ParamsDict["A_gnt"] = A_gnt
    ParamsDict["D_nt"] = D_nt
    ParamsDict["AH_nt"] = AH_nt
    ParamsDict["AR_nt"] = AR_nt

    if nosun
        gsun = 3
        ParamsDict["A_gnt"][gsun,:,:] .= 0
    end

    # Forming parameters struct
    parameters = Params(
        ParamsDict["region_n"], ParamsDict["max_dem_n"], ParamsDict["technology_g"], ParamsDict["G"], ParamsDict["G_r"], ParamsDict["N"], ParamsDict["L"], ParamsDict["L_ind"], 
        ParamsDict["T"], ParamsDict["S"], ParamsDict["H"], ParamsDict["κ"], ParamsDict["μ"], ParamsDict["C"], ParamsDict["C̄"], ParamsDict["C_E"], ParamsDict["R_E"], ParamsDict["τ_t"],
        ParamsDict["Gmin_gn"], ParamsDict["Gmax_gn"], ParamsDict["A_gnt"], ParamsDict["D_nt"], ParamsDict["I_g"], ParamsDict["M_g"], ParamsDict["C_g"], ParamsDict["e_g"], ParamsDict["E_g"],
        ParamsDict["r⁻_g"], ParamsDict["r⁺_g"], ParamsDict["I_l"], ParamsDict["M_l"], ParamsDict["C_l"], ParamsDict["B_l"], ParamsDict["e_l"], ParamsDict["Tmin_l"], ParamsDict["Tmax_l"],
        ParamsDict["ξ_s"], ParamsDict["I_s"], ParamsDict["C_s"], ParamsDict["Smin_sn"], ParamsDict["Smax_sn"], ParamsDict["Wmax_hn"], ParamsDict["Wmin_hn"], ParamsDict["Hmax_hn"],
        ParamsDict["Hmin_hn"], ParamsDict["HRmax_n"], ParamsDict["Fmin_n"], ParamsDict["AH_nt"], ParamsDict["AR_nt"], ParamsDict["I_h"], ParamsDict["M_h"], ParamsDict["C_h"],
        ParamsDict["e_h"], ParamsDict["E_h"], ParamsDict["r⁻_h"], ParamsDict["r⁺_h"]
    );
    return parameters
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
flatten(x::AffExpr) = x
flatten(x::VecOrMat{AffExpr}) = x[:]
flatten(x::Array{<:Any}) = x[:]

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