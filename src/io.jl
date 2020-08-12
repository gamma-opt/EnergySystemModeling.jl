using CSV, JSON, DataFrames, JLD, MAT

"""Equivalent annual cost (EAC)

# Arguments
* `cost::Real`: Net present value of a project.
* `n::Integer`: Number of payments.
* `r::Real`: Interest rate.
"""
function equivalent_annual_cost(cost::Real, n::Integer, r::Real)
    factor = if iszero(r) n else (1 - (1 + r)^(-n)) / r end
    cost / factor
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
- `instance_path::AbstractString`: Path to the instance directory.
"""
function Params(instance_path::AbstractString) 
    # Load indexes and constant parameters
    indices = JSON.parsefile(joinpath(instance_path, "indices.json"))

    # TODO: implement time period clustering: τ, T, τ_t

    # Load indices. Convert JSON values to right types.
    G = indices["G"] |> Array{Int}
    G_r = indices["G_r"] |> Array{Int}
    N = indices["N"] |> Array{Int}
    L = indices["L"] |> Array{Array{Int}}
    τ = 1
    T = 1:indices["T"]
    S = indices["S"] |> Array{Int}

    # Load constant parameters
    constants = JSON.parsefile(joinpath(instance_path, "constants.json"))
    κ = constants["kappa"]
    C = constants["C"]
    C̄ = constants["C_bar"]
    interest_rate = constants["r"]

    # Load time clustered parameters
    τ_t = ones(length(T))
    Q_gn = zeros(length(G), length(N))
    D_nt = zeros(length(N), length(T))
    A_gnt = ones(length(G), length(N), length(T))
    W_nmax = zeros(length(N))
    W_nmin = zeros(length(N))
    f_int = zeros(length(N), length(T))
    f′_int = zeros(length(N), length(T))
    H_n = zeros(length(N))
    H′_n = zeros(length(N))
    F_onmin = zeros(length(N))
    for n in N
        # Load node values from CSV files.
        df = CSV.read(joinpath(instance_path, "nodes", "$n.csv")) |> DataFrame
        capacitydf = CSV.read(joinpath(instance_path, "capacity.csv")) |> DataFrame
        for g in G
            Q_gn[g, n] = capacitydf[n, g+1]
        end
        D_nt[n, :] = df.Demand
        A_gnt[1, n, :] = df.Avail_Wind_On
        A_gnt[2, n, :] = df.Avail_Wind_Off
        A_gnt[3, n, :] = df.Avail_Sol
        A_gnt[A_gnt .< 0.01] .= 0
        W_nmax[n] = capacitydf.Max_Hyd_Level[n]
        W_nmin[n] = capacitydf.Min_Hyd_Level[n]
        f_int[n,:] = df.Hyd_In
        f′_int[n,:] = df.HydRoR_In
        H_n[n] = capacitydf.Hydro[n]
        H′_n[n] = capacitydf.HydroRoR[n]
        F_onmin[n] = (sum(f_int[n, :]) / 8760) * 0.05
        
    end

    # Load technology parameters
    technology = joinpath(instance_path, "technology.csv") |>
        CSV.read |> DataFrame
    I_g = equivalent_annual_cost.(technology.investment_cost .* 1000, technology.lifetime,
                                  interest_rate)
    M_g = technology.fixedOM .* 1000
    C_g = technology.fuel_cost ./ technology.efficiency .+ technology.varOM
    r⁻_g = technology.r_minus
    r⁺_g = technology.r_plus

    # Load transmission parameters
    transmission = joinpath(instance_path, "transmission.csv") |>
        CSV.read |> DataFrame
    M_l = transmission.M
    I_l = equivalent_annual_cost.(transmission.cost .* transmission.dist .+ M_l,
                                  transmission.lifetime, interest_rate)
    C_l = transmission.C
    B_l = transmission.B

    # Load storage parameters
    storage = joinpath(instance_path, "storage.csv") |>
        CSV.read |> DataFrame
    ξ_s = storage.xi
    I_s = equivalent_annual_cost.(storage.cost, storage.lifetime, interest_rate)
    C_s = storage.C
    b0_sn = storage[:, [Symbol("b0_$n") for n in N]] |> Matrix

    # Return Params struct
    Params(
        G, G_r, N, L, T, S, κ, C, C̄, τ, τ_t, Q_gn, A_gnt, D_nt, I_g, M_g, C_g,
        r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b0_sn,
        W_nmax, W_nmin, f_int, f′_int, H_n, H′_n, F_onmin)
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

function convert_type(::Type{Array{T, N}}) where T <: Number where N
    t = T
    for i in 1:N
        t = Array{t, 1}
    end
    return t
end
convert_type(::Type{Array{T, 1}}) where T <: Array = Array{T, 1}
convert_type(t::Type{T}) where T <: Number = t

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
function replace_nans(array::Array{Float64, N}) where N
    for i = eachindex(array)
        if isnan(array[i])
            array[i] = zero(1)
        end
    end
end

"""Reads wind, solar and hydro data produced by the GlobalEnergyGIS package into CSV node files.

# Arguments
- `instance_path::AbstractString`: Path to the instance directory.
- `era_year::AbstractString`: Year of the ERA5 data used in the GlobalEnergyGIS package.
- `era_year::AbstractString`: Name of the GIS region used.
"""
function create_nodedata(instance_path::AbstractString, era_year::AbstractString, gisregion::AbstractString)

    #Read files
    solarvars = matread(joinpath(instance_path, "GISdata_solar$(era_year)_$gisregion.mat"))
    windvars = matread(joinpath(instance_path, "GISdata_wind$(era_year)_$gisregion.mat"))
    hydrovars = matread(joinpath(instance_path, "GISdata_hydro_$gisregion.mat"))
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
    replace_nans(CFtime_cspplantA)
    replace_nans(CFtime_cspplantB)
    replace_nans(CFtime_pvplantA)
    replace_nans(CFtime_pvplantB)
    replace_nans(CFtime_pvrooftop)

    #Preallocate arrays for results
    avail_sol_cspA = zeros(Float64, 8760, n)
    avail_sol_cspB = zeros(Float64, 8760, n)
    avail_sol_pvA = zeros(Float64, 8760, n)
    avail_sol_pvB = zeros(Float64, 8760, n)
    avail_sol_pvrooftop = zeros(Float64, 8760, n)

    #Calculate absolute availability values and sum up the different classes for each solar power type
    for i in 1:8760
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
    replace_nans(CFtime_windoffshore)
    replace_nans(CFtime_windonshoreA)
    replace_nans(CFtime_windonshoreB)

    #Preallocate arrays for results
    avail_wind_offshore = zeros(Float64, 8760, n)
    avail_wind_onshoreA = zeros(Float64, 8760, n)
    avail_wind_onshoreB = zeros(Float64, 8760, n)

    #Calculate absolute availability values and sum up the different classes for each wind power type
    for i in 1:8760
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
    replace_nans(existinginflow)
    #Preallocate array for intermediate result
    avail_inflow = zeros(Float64, 8760, n)

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
        nodedata = zeros(8760, 6)
        nodedata[:,1] = demandvars[:,i]
        nodedata[:,2] = avail_sol[:,i]
        nodedata[:,3] = avail_wind_on[:,i]
        nodedata[:,4] = avail_wind_off[:,i]
        nodedata[:,5] = hyd_in[:,i]
        nodedata[:,6] = hydRoR_in[:,i]
        nodedata = convert(DataFrame, nodedata)
        rename!(nodedata, ["Demand", "Avail_Sol", "Avail_Wind_On", "Avail_Wind_Off", "Hyd_In", "HydRoR_In"])
        CSV.write(joinpath(instance_path, "nodes", "$i.csv"), nodedata)
    end

    
end

function getdispatch(output_path::AbstractString)

    variables = JSON.parsefile(joinpath(output_path, "variables.json"))
    parameters = JSON.parsefile(joinpath(output_path, "parameters.json"))

    pgnt = variables["p_gnt"] |> Array{Array{Array{Float64}}}
    hyd = variables["h_nt"] |> Array{Array{Float64}}
    dem = parameters["D_nt"] |> Array{Array{Float64}}
    dispatch = zeros(13, 11)
    for n in 1:11
        for g in 1:8
            for t in 1:8760
                dispatch[n, g] = dispatch[n, g] + pgnt[t][n][g]
            end  
        end
        for t in 1:8760
            dispatch[n, 9] = dispatch[n, 9] + hyd[t][n]
            dispatch[n, 11] = dispatch[n, 11] + dem[t][n]
        end
        dispatch[n, 10] = sum(dispatch[n, 1:9])
    end
    for i in 1:11
        dispatch[12, i] = sum(dispatch[:, i])
    end
    dispatch[13, :] = dispatch[12, :] ./ dispatch[12, 10]
    #
    #capacity = variables["p̄_gn"] |> Array{Array{Float64}}
    #hydcap = variables["H_n"] |> Array{Float64}
    #hydRoRcap = variables["H′_n"] |> Array{Float64}
    #capacities = zeros(13, 10)
    #for n in 1:11
    #    for g in 1:8
    #        capacities[n, g] = capacity[n, g]
    #   end
    #    for t in 1:8760
    #        dispatch[n, 9] = dispatch[n, 9] + hyd[t][n]
    #    end
    #    dispatch[n, 10] = sum(dispatch[n, :])
    #end
    #for i in 1:10
    #    dispatch[12, i] = sum(dispatch[:, i])
    #end
    #dispatch[13, :] = dispatch[12, :] ./ dispatch[12, 10]

    dispatch = convert(DataFrame, dispatch)
    rename!(dispatch, ["WIND_ON", "WIND_OFF", "SOLAR", "BIOMASS", "NUCLEAR", "COAL", "GAS_CC", "GAS_OC", "HYDRO", "TOTAL", "DEMAND"])
    CSV.write(joinpath(output_path, "dispatch.csv"), dispatch)

end