
Europe = [
    # "Country"                   "GADM"    
      "Albania"                   GADM("Albania") 
      "Armenia"                   GADM("Armenia") 
      "Austria"                   GADM("Austria") 
      "Azerbaijan"                GADM("Azerbaijan") 
      "Belarus"                   GADM("Belarus")                
      "Belgium"                   GADM("Belgium") 
      "Bosnia and Herzegovina"    GADM("Bosnia and Herzegovina")            
      "Bulgaria"                  GADM("Bulgaria")           
      "Croatia"                   GADM("Croatia")        
      "Cyprus"                    GADM("Cyprus")             
      "Czech Republic"            GADM("Czech Republic")     
      "Denmark"                   GADM("Denmark")          
      "Estonia"                   GADM("Estonia")           
      "Finland"                   GADM("Finland")         
      "France"                    GADM("France")              
      "Germany"                   GADM("Germany")         
      "Greece"                    GADM("Greece")           
      "Hungary"                   GADM("Hungary")     
      "Iceland"                   GADM("Iceland")             
      "Ireland"                   GADM("Ireland")           
      "Italy"                     GADM("Italy")  
      "Kosovo"                    GADM("Kosovo")                        
      "Latvia"                    GADM("Latvia")           
      "Lithuania"                 GADM("Lithuania")        
      "Luxembourg"                GADM("Luxembourg")        
      "Malta"                     GADM("Malta")     
      "Moldova"                   GADM("Moldova") 
      "Montenegro"                GADM("Montenegro")                              
      "Netherlands"               GADM("Netherlands")   
      "North Macedonia"           GADM("Macedonia")    
      "Norway"                    GADM("Norway")          
      "Poland"                    GADM("Poland")             
      "Portugal"                  GADM("Portugal")           
      "Romania"                   GADM("Romania") 
      "Russia"                    GADM("Russia")
      "Serbia"                    GADM("Serbia")              
      "Slovakia"                  GADM("Slovakia")           
      "Slovenia"                  GADM("Slovenia")          
      "Spain"                     GADM("Spain")              
      "Sweden"                    GADM("Sweden")            
      "Switzerland"               GADM("Switzerland")  
      "Turkey"                    GADM("Turkey")
      "Ukraine"                   GADM("Ukraine")         
      "United Kingdom"            GADM("United Kingdom")     
  ]  



Regions_dict = Dict( "Nordics" => (["Finland", "Sweden", "Norway", "Denmark"]),
  "Central" => (["Germany", "Austria", "Switzerland", "Czech Republic"]),
  "Western" => (["France", "United Kingdom", "Ireland", "Netherlands", "Belgium", "Luxembourg"]),
  "Mediterranian" => (["Spain", "Portugal", "Italy", "Greece", "Croatia", "Malta", "Albania", "Bosnia and Herzegovina"]),
  "Eastern" => (["Poland", "Slovakia", "Hungary", "Lithuania", "Latvia", "Estonia"]))


"""Takes the list of countries specified in run_data_generation.jl and turns them into a format that GlobalEnergyGIS can work with.
# Arguments
- `Regions`: List of regions used for the model. 
"""
# function if one wants specific countries
# Check if the regions specified in run_data_generation.jl are in Europe or in Regions_dict
# Get corresponding values
function get_countries(Regions)
    if sum(occursin.(Regions[1,1], Europe[:,1]))>0
        n = length(Europe[:,1])
        m = length(Regions)
        Dataset_Countries = []
        GADM_Name = []
        for j in 1:m
            for i in 1:n
                    if occursin.(Regions[j,1], Europe[i,1])
                        Dataset_Countries = [Dataset_Countries; permutedims(Regions[j,:])]
                        GADM_Name = [GADM_Name; GADM(Europe[i])]
                    end
            end
        end
        Countries = hcat(Dataset_Countries, GADM_Name)
    elseif sum(occursin.(Regions[1,1], collect(keys(Regions_dict))))>0
        Keys = []
        Values = []
        m = length(Regions)
        R = getindex.(Ref(Regions_dict),(Regions))
        for i in 1:m
            RVal = getindex.(Ref(Regions_dict),(Regions))[i]
            Keys = [Keys; Regions[i]]
            Values = [Values; GADM(RVal...)]
        end
        GADM_List = hcat(Keys, Values)
    end
end


# function if one wants to create whole regions
# if you want all countries in Europe
function get_countries()
   Europe
end


function replace_nans!(array::Array{Float64, N}) where N
    for i in eachindex(array)
        if isnan(array[i])
            array[i] = zero(1)
        end
    end
end

"""Uses wind, solar and hydro data produced by the GlobalEnergyGIS package to generate CSV files with the data used in the model.
# Arguments
- `inputdata`: Path to the output folder of GlobalEnergyGIS.
- `sspscenario_input`: Specified ssp scenario used to generate data.
- `sspyear_input`: Target year of the model.
- `era_year_input`: Year of the base data.
- `Dataset`: Countries generated by the function `get_countries`.
- `folder`: Sets the folder in which the generated data is put.
- `subfolder`: Sets the subfolder in which the generated data is put.
- `instance`: Sets the name of the instance for which the data is being generated.
- `T`: Time period of the data set in hours.
- `t`: number of technologies used in the model
"""

function create_data_sets(inputdata, sspscenario_input, sspyear_input, era_year_input, Dataset, folder, subfolder, instance, T, t, Fmin)

    # create path for the instance
    structure_path = joinpath(folder,subfolder)
    instance_path = mkpath(joinpath(structure_path, instance))
    nodes_path =  mkpath(joinpath(instance_path, "nodes"))

# Generate datasets for regions using GlobalEnergyGIS
    # define regions and countries and name the instance
 
   saveregions(instance, Dataset)
   makedistances(instance)
   createmaps(instance)

    # generate VRE data and demand
   GISsolar(gisregion = instance)
   GISwind(gisregion = instance)
   GIShydro(gisregion = instance)
   predictdemand(gisregion = instance, sspscenario=sspscenario_input, sspyear=sspyear_input, era_year=era_year_input)

# function get data
    # Name of the instance you defined in Inputdata.jl
    SolarData = matread(joinpath(inputdata, "GISdata_solar$(era_year_input)_$instance.mat"))
    WindData = matread(joinpath(inputdata, "GISdata_wind$(era_year_input)_$instance.mat"))
    HydroData = matread(joinpath(inputdata, "GISdata_hydro_$instance.mat"))
    TransmissionData = matread(joinpath(inputdata, "distances_$instance.mat"))

    Regionlist = typeof(TransmissionData["regionlist"]) == Float64 ? [TransmissionData["regionlist"]] : TransmissionData["regionlist"]
    n = length(Regionlist)

    # Solar

    PVCapA = typeof(SolarData["capacity_pvplantA"]) == Float64 ? [SolarData["capacity_pvplantA"]] : SolarData["capacity_pvplantA"]
    PVCapB = typeof(SolarData["capacity_pvplantB"]) == Float64 ? [SolarData["capacity_pvplantB"]] : SolarData["capacity_pvplantB"]
    CSPCapA = typeof(SolarData["capacity_cspplantA"]) == Float64 ? [SolarData["capacity_cspplantA"]] : SolarData["capacity_cspplantA"]
    CSPCapB = typeof(SolarData["capacity_cspplantB"]) == Float64 ? [SolarData["capacity_cspplantB"]] : SolarData["capacity_cspplantB"]
    RooftopCap = typeof(SolarData["capacity_pvrooftop"]) == Float64 ? [SolarData["capacity_pvrooftop"]] : SolarData["capacity_pvrooftop"]
    SolarCapTotal = PVCapA .+ PVCapB .+ CSPCapA .+ CSPCapB .+ RooftopCap
    replace_nans!(PVCapA)
    replace_nans!(PVCapB)
    replace_nans!(CSPCapA)
    replace_nans!(CSPCapB)
    replace_nans!(RooftopCap)
    replace_nans!(SolarCapTotal)

    CFtime_PVA = typeof(SolarData["CFtime_pvplantA"]) == Float64 ? [SolarData["CFtime_pvplantA"]] : SolarData["CFtime_pvplantA"]
    CFtime_PVB = typeof(SolarData["CFtime_pvplantB"]) == Float64 ? [SolarData["CFtime_pvplantB"]] : SolarData["CFtime_pvplantB"]
    CFtime_CSPA = typeof(SolarData["CFtime_cspplantA"]) == Float64 ? [SolarData["CFtime_cspplantA"]] : SolarData["CFtime_cspplantA"]
    CFtime_CSPB = typeof(SolarData["CFtime_cspplantB"]) == Float64 ? [SolarData["CFtime_cspplantB"]] : SolarData["CFtime_cspplantB"]
    CFtime_PVrooftop = typeof(SolarData["CFtime_pvrooftop"]) == Float64 ? [SolarData["CFtime_pvrooftop"]] : SolarData["CFtime_pvrooftop"]
    CFtime_solar = CFtime_PVA .+ CFtime_PVB .+ CFtime_CSPA .+ CFtime_CSPB .+ CFtime_PVrooftop
    replace_nans!(CFtime_PVA)
    replace_nans!(CFtime_PVB)
    replace_nans!(CFtime_CSPA)
    replace_nans!(CFtime_CSPB)
    replace_nans!(CFtime_PVrooftop)
    replace_nans!(CFtime_solar)

    #Preallocate arrays for results
    avail_sol_cspA = zeros(Float64, T, n)
    avail_sol_cspB = zeros(Float64, T, n)
    avail_sol_pvA = zeros(Float64, T, n)
    avail_sol_pvB = zeros(Float64, T, n)
    avail_sol_pvrooftop = zeros(Float64, T, n)

    #Calculate absolute availability values and sum up the different classes for each solar power type
    for i in 1:T
        avail_sol_cspA[i,:,:] = permutedims(sum(CSPCapA .* CFtime_CSPA[i,:,:], dims=2))
        avail_sol_cspB[i,:,:] = permutedims(sum(CSPCapB .* CFtime_CSPB[i,:,:], dims=2))
        avail_sol_pvA[i,:,:] = permutedims(sum(PVCapA .* CFtime_PVA[i,:,:], dims=2))
        avail_sol_pvB[i,:,:] = permutedims(sum(PVCapB .* CFtime_PVB[i,:,:], dims=2))
        avail_sol_pvrooftop[i,:,:] = permutedims(sum(RooftopCap .* CFtime_PVrooftop[i,:,:], dims=2))
    end
    #Calculate total solar capacity
    capacity_solar = permutedims(sum(CSPCapA, dims=2)) + permutedims(sum(CSPCapB, dims=2)) + permutedims(sum(PVCapA, dims=2)) +
                    permutedims(sum(PVCapB, dims=2))  + permutedims(sum(RooftopCap, dims=2))
    #Sum the different solar types and divide with total capacity to get relative availability values
    avail_sol = (avail_sol_cspA + avail_sol_cspB + avail_sol_pvA + avail_sol_pvB + avail_sol_pvrooftop) ./ capacity_solar

    # Wind
    WindCapA = typeof(WindData["capacity_onshoreA"]) == Float64 ? [WindData["capacity_onshoreA"]] : WindData["capacity_onshoreA"]
    WindCapB = typeof(WindData["capacity_onshoreB"]) == Float64 ? [WindData["capacity_onshoreB"]] : WindData["capacity_onshoreB"]
    WindCapTotal = WindCapA .+ WindCapB
    replace_nans!(WindCapA)
    replace_nans!(WindCapB)
    replace_nans!(WindCapTotal)

    WindCapOff = typeof(WindData["capacity_offshore"]) == Float64 ? [WindData["capacity_offshore"]] : WindData["capacity_offshore"]
    replace_nans!(WindCapOff)

    ## Wind Onshore
    CFtime_windonshoreA = typeof(WindData["CFtime_windonshoreA"]) == Float64 ? [WindData["CFtime_windonshoreA"]] : WindData["CFtime_windonshoreA"]
    CFtime_windonshoreB = typeof(WindData["CFtime_windonshoreB"]) == Float64 ? [WindData["CFtime_windonshoreB"]] : WindData["CFtime_windonshoreB"]
    CFtime_windonshore = CFtime_windonshoreA .+ CFtime_windonshoreB
    replace_nans!(CFtime_windonshoreA)
    replace_nans!(CFtime_windonshoreB)
    replace_nans!(CFtime_windonshore)

    ## Wind Offshore
    CFtime_windoffshore = typeof(WindData["CFtime_windoffshore"]) == Float64 ? [WindData["CFtime_windoffshore"]] : WindData["CFtime_windoffshore"]
    replace_nans!(CFtime_windoffshore)
 
    #Preallocate arrays for results
    avail_wind_offshore = zeros(Float64, T, n)
    avail_wind_onshoreA = zeros(Float64, T, n)
    avail_wind_onshoreB = zeros(Float64, T, n)

    # Calculate absolute availability values and sum up the different classes for each wind power type
    for i in 1:T
        avail_wind_offshore[i,:,:] = permutedims(sum(WindCapOff .* CFtime_windoffshore[i,:,:], dims=2))
        avail_wind_onshoreA[i,:,:] = permutedims(sum(WindCapA .* CFtime_windonshoreA[i,:,:], dims=2))
        avail_wind_onshoreB[i,:,:] = permutedims(sum(WindCapB .* CFtime_windonshoreB[i,:,:], dims=2))
    end
    capacity_wind_on = permutedims(sum(WindCapA, dims=2)) + permutedims(sum(WindCapB, dims=2))
    capacity_wind_off = permutedims(sum(WindCapOff, dims=2))
    avail_wind_on = (avail_wind_onshoreA + avail_wind_onshoreB) ./ capacity_wind_on
    avail_wind_off = avail_wind_offshore ./ capacity_wind_off
    replace_nans!(avail_wind_off)
    
      
    # demand
    gisdemand = JLD.load(joinpath(inputdata, "SyntheticDemand_$(instance)_$(sspscenario_input)-$(sspyear_input)_$(era_year_input).jld"), "demand")

    demand = zeros(Float64, T, n)
    for i in 1:T
        demand[i,:] = gisdemand[i,:]
    end

    # Transmission distances and connected nodes
    Distances = typeof(TransmissionData["distances"]) == Float64 ? [TransmissionData["distances"]] : TransmissionData["distances"]
    connected = typeof(TransmissionData["connected"]) == Float64 ? [TransmissionData["connected"]] : TransmissionData["connected"]

    # Hydro
    HydroCap = typeof(HydroData["existingcapac"]) == Float64 ? [HydroData["existingcapac"]] : HydroData["existingcapac"]
    replace_nans!(HydroCap)
   
    existingcapac = permutedims(HydroData["existingcapac"]) .* 1000
    potentialcapacities = HydroData["potentialcapac"]
    potentialcapac = []
    for i in 1:n
        potentialcapac = vcat(potentialcapac, sum(potentialcapacities[i,:,:])*1000 + existingcapac[1,i])
    end
    #Monthly inflow
    existinginflow = permutedims(HydroData["existinginflowcf"])
    # Replace NaNs from monthly data with 0's
    replace_nans!(existinginflow)
    # Preallocate array for intermediate result
    avail_inflow = zeros(Float64, T, n)           ## Make it variable

    #Turn monthly inflow data into hourly data
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    lasthours = 24 * cumsum(days)
    firsthours = [1; 1 .+ lasthours[1:end-1]]

    for m in 1:12
        if lasthours[m] <= T
            for i = firsthours[m]:lasthours[m]
                avail_inflow[i,:] = existinginflow[m,:]
            end
        elseif lasthours[m] >= T
            for i = firsthours[m]:T
                avail_inflow[i,:] = existinginflow[m,:]
            end
        end
    end
      
    # Get the installed capacity of Run of river and Hydro reservoir plants from ENTSO-E datasets: https://transparency.entsoe.eu/generation/r2/installedGenerationCapacityAggregation/show
    # Create a table (RoR_Res) with values from ENTSOE for Pumped Storage, Run of River and Reservoir for each values
    data_path = joinpath("ENTSO-E_data", "Capacities")
    if sum(occursin.(Regions[1,1], Europe[:,1]))>0
        Countries = Regions[:,1]
        Hydro_Matrix = Array{Any}(undef, 0, 3)
        # take the values from the ENTSO-E CSV files, put them into a matrix
        for i in 1:length(Countries)
            Country = Regions[:,1][i]
            Country_ENTSOE = joinpath(data_path, "$(Country).csv") |> CSV.File |> DataFrame
            RoR_Res_data = Country_ENTSOE[in(["Hydro Pumped Storage", "Hydro Run-of-river and poundage", "Hydro Water Reservoir"]).(Country_ENTSOE."Production Type"), 2]
            replace!(RoR_Res_data, "n/e"=>"0")
            RoR_Res_data = parse.(Float64, string.(RoR_Res_data))
            Water = RoR_Res_data
            Hydro_Matrix = vcat(Hydro_Matrix, permutedims(Water))
        end
        Hydro = hcat(Countries, Hydro_Matrix)
        Header = ["Country" "Pumped Storage" "Run of River" "Reservoir"]
        RoR_Res = vcat(Header, Hydro)
        replace!(RoR_Res, "N/A"=>0)
        # Reservoir percentage is the reservoir value devided by the sum of all three hydro power generation capacities
        reservoir = RoR_Res[2:end,4]./(RoR_Res[2:end,2] .+ RoR_Res[2:end,3] .+ RoR_Res[2:end,4])
        replace_nans!(reservoir)
        reservoir_table = hcat(Regions, reservoir)
        ## Percentage between PHS and reservoir for later calculation of Wmin and Wmax
        reservoir_PHS = RoR_Res[2:end,4]./(RoR_Res[2:end,2] .+ RoR_Res[2:end,4])
        replace_nans!(reservoir_PHS)
    elseif sum(occursin.(Regions[1,1], collect(keys(Regions_dict))))>0
        Hydro_Matrix = Array{Any}(undef, 0, 3)
        ## get the countries, that are in the regions
        Countries = []
        for i in 1:length(Regions)
            R = getindex.(Ref(Regions_dict),(Regions))[i]
            Countries = [Countries; R]
        end
        Hydro_Matrix = Array{Any}(undef, 0, 3)
        ## get ENTSO-E values for the countries in the regions
        for i in 1:length(Countries)
            Country = Countries[:,1][i]
            Country_ENTSOE = joinpath(data_path, "$(Country).csv") |> CSV.File |> DataFrame
            RoR_Res_data = Country_ENTSOE[in(["Hydro Pumped Storage", "Hydro Run-of-river and poundage", "Hydro Water Reservoir"]).(Country_ENTSOE."Production Type"), 2]
            replace!(RoR_Res_data, "n/e"=>"0")
            replace!(RoR_Res_data, "N/A"=>"0")
            RoR_Res_data = parse.(Float64, string.(RoR_Res_data))
            Water = RoR_Res_data
            Hydro_Matrix = vcat(Hydro_Matrix, permutedims(Water))
        end
        Hydro = hcat(Countries, Hydro_Matrix)
        Header = ["Country" "Pumped Storage" "Run of River" "Reservoir"]
        RoR_Res = vcat(Header, Hydro)
        replace!(RoR_Res, "n/e"=>0)
        replace!(RoR_Res, "N/A"=>0)
        RoR_Res[2:end,2:end] = parse.(Float64, string.(RoR_Res[2:end,2:end]))   
        res_table = []
        ## get the region of each country
        for i in 1:length(Regions)
            R = getindex.(Ref(Regions_dict),(Regions))[i]
            Key = [k for (k,v) in Regions_dict if v==R]
            Key_List = repeat(Key, inner=length(R))
            res_table = [res_table; Key_List]
        end
        # List all countries with their region
        List_full = DataFrame(hcat(res_table, RoR_Res[2:end,2:end]), :auto)
        # sum the different regions, so you get one value per region
        Hydro_summed = combine(groupby(List_full, :x1), names(List_full, Not(:x1)) .=> sum, renamecols=false)
        reservoir = Hydro_summed[!,:x4]./(Hydro_summed[!,:x2] .+ Hydro_summed[!,:x3] .+ Hydro_summed[!,:x4])
        reservoir_table = hcat(Regions, reservoir)        
        ## Percentage between PHS and reservoir for later calculation of Wmin and Wmax
        reservoir_PHS = Hydro_summed[!,:x4]./(Hydro_summed[!,:x2] .+ Hydro_summed[!,:x4])
        replace_nans!(reservoir_PHS)
    end   

    reservoirp =  permutedims(reservoir_table[:,2])

    #Hydro capacities
    hydroCap = existingcapac .* reservoirp
    hydroRoRCap = existingcapac - hydroCap
    hydrocapacity = [hydroCap; hydroRoRCap]
    
    #Reservoir and RoR inflow
    hyd_in = avail_inflow .* hydroCap           ## availability is not what is installed, but the amount of water available
    hydRoR_in = avail_inflow .* hydroRoRCap


# Hydro_capacity.csv file

    ## get the reservoir levels from ENTSO-E, base year 2018: https://transparency.entsoe.eu/generation/r2/waterReservoirsAndHydroStoragePlants/show
    ## Wmax is the maximum value of year 2018 for each country
    ## Wmin is the minimum value of year 2018 for each country
    data_path_reslvl = joinpath("ENTSO-E_data", "Reservoir_level")
    if sum(occursin.(Regions[1,1], Europe[:,1]))>0
        Countries = Regions[:,1]
        Hydro_Matrix = Array{Any}(undef, 0, 3)
        W_Min = Array{Any}(undef, 0, 1)
        W_Max = Array{Any}(undef, 0, 1)
        # take the values from the ENTSO-E CSV files, put them into a matrix
        for i in 1:length(Countries)
            Country = Regions[:,1][i]
            Country_ENTSOE = joinpath(data_path_reslvl, "$(Country).csv") |> CSV.File |> DataFrame
            W_data_total = Matrix(Country_ENTSOE[!, Not("Week")])
            replace!(W_data_total, "n/e"=>"0")
            replace!(W_data_total, "N/A"=>"0")
            W_data_total = parse.(Float64, string.(W_data_total))
            W_Min = vcat(W_Min, minimum(W_data_total, dims=1))
            W_Max = vcat(W_Max, maximum(W_data_total, dims=1))
        end   
        Wmin_summed = W_min_total = hcat(Countries, W_Min)
        Wmax_summed = W_max_total = hcat(Countries, W_Max) 
    elseif sum(occursin.(Regions[1,1], collect(keys(Regions_dict))))>0
        ## get the countries, that are in the regions
        Countries = []
        for i in 1:length(Regions)
            R = getindex.(Ref(Regions_dict),(Regions))[i]
            Countries = [Countries; R]
        end
        W_Min = Array{Any}(undef, 0, 1)
        W_Max = Array{Any}(undef, 0, 1)
        for i in 1:length(Countries)
            Country = Countries[:,1][i]
            Country_ENTSOE = joinpath(data_path_reslvl, "$(Country).csv") |> CSV.File |> DataFrame
            W_data_total = Matrix(Country_ENTSOE[!, Not("Week")])
            replace!(W_data_total, "n/e"=>"0")
            replace!(W_data_total, "N/A"=>"0")
            W_data_total = parse.(Float64, string.(W_data_total))
            W_Min = vcat(W_Min, minimum(W_data_total, dims=1))
            W_Max = vcat(W_Max, maximum(W_data_total, dims=1))
        end
        W_min_total = hcat(Countries, W_Min)
        W_max_total = hcat(Countries, W_Max)
        W_table = []
        for i in 1:length(Regions)
            R = getindex.(Ref(Regions_dict),(Regions))[i]
            Key = [k for (k,v) in Regions_dict if v==R]
            Key_List = repeat(Key, inner=length(R))
            W_table = [W_table; Key_List]
        end
        Wmin_full = DataFrame(hcat(W_table, W_min_total[:,2]), :auto)
        Wmin_summed = combine(groupby(Wmin_full, :x1), names(Wmin_full, Not(:x1)) .=> sum, renamecols=false)
        Wmax_full = DataFrame(hcat(W_table, W_max_total[:,2]), :auto)
        Wmax_summed = combine(groupby(Wmax_full, :x1), names(Wmax_full, Not(:x1)) .=> sum, renamecols=false)
    end    

    H_max = potentialcapac .* permutedims(reservoirp) # reservoir production capacity per node [MWh]
    NodeNr = [1:n;]
    H_min = permutedims(hydroCap)      # existing reservoir capacity per node (percentage)
    HR_max = permutedims(hydroRoRCap)        # RoR counterpart of existing capacity 
    ## ENTSO-E data includes both PHS and reservoir, so we need the reservoir percentage to calculate acutal Wmin and Wmax for reservoir only
    W_max = Wmax_summed[:,2].*reservoir_PHS
    W_min = Wmin_summed[:,2].*reservoir_PHS
    hydro_tech = repeat(1:1, inner=n)
    Hydro_cap = DataFrame(hcat(NodeNr, hydro_tech, H_min, H_max, W_min, W_max), :auto)
    rename!(Hydro_cap, ["node", "hydro_tech", "hcap_min", "hcap_max", "wcap_min" ,"wcap_max"])
    CSV.write(joinpath(instance_path, "hydro_capacity.csv"), Hydro_cap)

# hydro.csv file
    node = [1:n;]
    HydroRoR = permutedims(hydroRoRCap)
    Fmin_mat = fill(Fmin, n)
    hyd_in_avg = []
    for i in 1:n
        hyd_in_avg = vcat(hyd_in_avg, mean(hyd_in[:,i]))
    end
    hydror_avg = []
    for i in 1:n
        hydror_avg = vcat(hydror_avg, mean(hydRoR_in[:,i]))
    end
    hyd_flow_min = []
    for i in 1:n
       hyd_flow_min = vcat(hyd_flow_min, argmin((avail_inflow .* existingcapac)[:,i]))
    end
    Hydro = DataFrame(hcat(node, HydroRoR, Fmin_mat, hyd_in_avg, hydror_avg, hyd_flow_min), :auto)
    rename!(Hydro, ["node", "HydroRoR", "Fmin", "hyd_in_avg", "hydror_avg" ,"hyd_flow_min"])
    CSV.write(joinpath(instance_path, "hydro.csv"), Hydro)    

# function create_node_data()
    # Create a Vector for every technology spanning all nodes (at 5 nodes: 43.800 columns)
    Demand = []
    WindOff = []
    WindOn = []
    Solar = []
    hydro_in = []
    hydroRoR_in = []
    
    for i in 1:n
        Demand$i = sum(demand[:,i], dims = 2)
        append!(Demand, Demand$i)
    end
    for i in 1:n
        WindOff$i = sum(avail_wind_off[:,i], dims = 2)
        append!(WindOff, WindOff$i)
    end
    for i in 1:n
        WindOn$i = sum(avail_wind_on[:,i], dims = 2)
        append!(WindOn, WindOn$i)
    end
    for i in 1:n
        Solar$i = sum(avail_sol[:,i], dims = 2)
        append!(Solar, Solar$i)
    end
    for i in 1:n
        hyd_in$i = sum(hyd_in[:,i], dims = 2)
        append!(hydro_in, hyd_in$i)
    end
    for i in 1:n
        hydRoR_in$i = sum(hydRoR_in[:,i], dims = 2)
        append!(hydroRoR_in, hydRoR_in$i)
    end
    All = hcat(Demand, Solar, WindOn, WindOff, hydro_in, hydroRoR_in) 
    
    # Generate CSV files for every single node
    x = 1
    for i = 1:n
        NodeData = All[x:T*i,:]
        x = x+T
        N_m = DataFrame(NodeData, :auto)
        rename!(N_m, ["Demand", "Avail_Sol", "Avail_Wind_On" ,"Avail_Wind_Off", "Hyd_In" ,"HydRoR_In"])
        CSV.write(joinpath(nodes_path, "$i.csv"), N_m)
    end

    x = 1
    # create file for node_specs
    NodeDemand = []
    for i = 1:n
        Max_Node_Demand = maximum(Demand[x:T*i], dims=1)
        x = x+T 
        NodeDemand = append!(NodeDemand, Max_Node_Demand)
    end

    NodeNr = [1:n;]
    Max_Demand = NodeDemand
    nodes_specs = DataFrame(hcat(NodeNr, Regionlist, Max_Demand), :auto)
    rename!(nodes_specs, ["Node", "Name", "Max_Demand" ])
    CSV.write(joinpath(instance_path, "nodes_specs.csv"), nodes_specs)

# function gen_capacity() 
    # capacity of the remaining non-renewable energy sources fixed at 1000 GW (biomass, nuclear, coal, gas_cc, gas_oc)
    Capacity_else = 1000*1000
    Capacity = zeros(5,1)
    gcap_max = []
            
    # create a vector with the capacity of each energy source per node
    for i in 1:n
        Capacity$i = [sum(WindCapTotal[i,:]), sum(WindCapOff[i,:]), sum(SolarCapTotal[i,:]),  Capacity_else, Capacity_else, Capacity_else, Capacity_else, Capacity_else].*1000
        append!(gcap_max, Capacity$i)
    end
     
    # Zero matrix with 4 columns
    gen_capacity = zeros(n*t, 4)
    
    # Column with node numbers from 1 to node n, which repeats itself as often as we have technologies
    node = repeat(1:n, inner=t)
    # repeating the technologies 1 - 8 for every node
    gen_tech = repeat([1:t;], n)
    # fill the whole column of gcap_min with 0
    gcap_min = fill(0, (n*t,1))  
    Gen_capac = hcat(node, gen_tech, gcap_min, gcap_max)
    gen_capacity = DataFrame(Gen_capac, :auto)                                                  
    rename!(gen_capacity, ["node", "gen_tech", "gcap_min" ,"gcap_max"])
    CSV.write(joinpath(instance_path, "gen_capacity.csv"), gen_capacity)

# sto_capacity()
    st = 1          
    Max_capacity = 1000000000
    s = repeat([1:st;], n)
    node = repeat(1:n, inner=st)
    scap_min = zeros(n)
    scap_max = fill(Max_capacity, (n*st, 1))
    sto_cap = DataFrame(hcat(s, node, scap_min, scap_max), :auto)
    rename!(sto_cap, ["s", "node", "scap_min", "scap_max"])
    CSV.write(joinpath(instance_path, "sto_capacity.csv"), sto_cap)

# transmission()
    Dist = Distances[tril!(trues(size(Distances)), -1)]
    # Connected nodes
    Con = connected[tril!(trues(size(connected)), -1)]
    # Find distances for connected nodes only
    Dist_Con = hcat(Dist, Con)
    Dist_Con_df = DataFrame(Dist_Con, :auto)
    Dist_con = Dist_Con_df[Dist_Con_df[!,:x2].!=0,:]
    Distance_Connected = Matrix(Dist_con)[:,1]
    # Ged rid of double nodes by only using lower triangular matrix
    LinesLT = LowerTriangular(connected)
    # Find connections by getting the CartesianIndex of the connected nodes
    LinesCI = findall(x->x==1, LinesLT)
    Lines = hcat(getindex.(LinesCI, 1), getindex.(LinesCI, 2))
    LinesDoubledf = DataFrame(Lines, :auto)
    # Delete entries of connections between only one node
    Linesdf = LinesDoubledf[LinesDoubledf.x1 .!= LinesDoubledf.x2, :]
    LinesMatrix = Matrix(Linesdf)
    # Reshape into a 1d array
    LinesSingle = collect(eachrow(LinesMatrix))
    L_ind = [1:size(LinesMatrix)[1];]
      
    # Create Matrix for remaining parameters of transmission
    TransData = [
       # :cost       :converter_cost         :M      :C      :B      :efficiency         :lifetime       :tcap_min       :tcap_max
        400          150000                  0.02    0       1       0.95                40              0               1000000000
    ]
    
    # Fill the matrix according to the number of pairs in the system
    TransDatafull = repeat(TransData; outer=[size(LinesMatrix)[1]])
    # combine the tables, add header and create CSV file
    TransmissionData = DataFrame(hcat(L_ind, LinesSingle, Distance_Connected, TransDatafull), :auto)
    rename!(TransmissionData, ["L_ind" ,"line" ,"dist" , "cost", "converter_cost", "M" ,"C" , "B", "efficiency", "lifetime" ,"tcap_min" ,"tcap_max"])
    CSV.write(joinpath(instance_path, "transmission.csv"), TransmissionData)

# function generate_tech_tables()
    techtable = [
        :name          :g       :investment_cost    :fixedOM    :varOM      :fuel_cost      :efficiency     :emissions      :lifetime       :r_minus    :r_plus
        :wind_on        1       1473                35          0           0               1               0               25              1           1
        :wind_off       2       3800                80          0           0               1               0               25              1           1
        :solar          3       995                 26          0           0               1               0               25              1           1
        :biomass        4       2141                100         0           7               0.448           0.39            30              1           1
        :nuclear        5       5000                150         3           3               0.34            0               60              0.05        0.05
        :coal           6       1800                25          6           8               0.466           0.34            40              0.15        0.15
        :gas_cc         7       800                 20          4           22              0.615           0.20            30              0.3         0.3
        :gas_oc         8       400                 15          3           22              0.395           0.20            30              1           1
    ]
    CSV.write(joinpath(instance_path, "gen_technology.csv"), Tables.table(techtable), writeheader=false)
        
    ## sources investment cost
    # Wind onshore, Wind offshore, Solar, Biomass: https://www.irena.org/publications/2020/Jun/Renewable-Power-Costs-in-2019 
    # Nuclear, coal, gas_cc: https://www.econstor.eu/bitstream/10419/80348/1/757528015.pdf (rough average of the different estimations)
    # 



    # Hydro
    hydro_tech = [
        :name          :hydro_tech       :investment_cost    :fixedOM    :varOM      :fuel_cost      :efficiency     :emissions      :lifetime       :r_minus    :r_plus
        :hydropower     1                 1704                200          0           0               1               0               60              1           1
    ]
    CSV.write(joinpath(instance_path, "hydro_technology.csv"), Tables.table(hydro_tech), writeheader=false)
    ## sources investment cost
    # https://www.irena.org/publications/2020/Jun/Renewable-Power-Costs-in-2019 

        
    # Storage
    storage = [
        :s          :xi      :C    :cost    :lifetime        
        1           0.86     15     185     20
    ]
    CSV.write(joinpath(instance_path, "storage.csv"), Tables.table(storage), writeheader=false)

    # create weights.csv if hours = 8760
    if T == 8760
        weights = repeat([1:1;], 8760)
        df = DataFrame(col1 = weights)
        rename!(df, ["Weights"])
        CSV.write(joinpath(instance_path, "weights.csv"), df)
    end
   
    # indices.json
    G = [1:t;]
    G_r = [1:4;]
    N = [1:n;]
    L = LinesSingle
    L_ind = L_ind
    S = [1]
    H = [1]
    T = T
    indices = Dict("G" => G, "G_r" => G_r, "N" => N, "L" => L, "L_ind" => L_ind, "S" => S, "H" => H, "T" => T)
    open(joinpath(instance_path, "indices.json"), "w") do f
        JSON.print(f, indices)
    end
 
end

