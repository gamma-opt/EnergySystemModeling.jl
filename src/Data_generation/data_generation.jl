function replace_nans!(array::Array{Float64, N}) where N
    for i in eachindex(array)
        if isnan(array[i])
            array[i] = zero(1)
        end
    end
end


function get_data(inputdata, regionset, sspscenario, sspyear, era_year)

    # Name of the regionset you defined in Inputdata.jl
    SolarData = matread(joinpath(inputdata, "GISdata_solar2018_$regionset.mat"))
    WindData = matread(joinpath(inputdata, "GISdata_wind2018_$regionset.mat"))
    HydroData = matread(joinpath(inputdata, "GISdata_hydro_$regionset.mat"))
    TransmissionData = matread(joinpath(inputdata, "distances_$regionset.mat"))

    # Get the different capacities per VRE source from GlobalEnergyGIS files
    WindCapA = typeof(WindData["capacity_onshoreA"]) == Float64 ? [WindData["capacity_onshoreA"]] : WindData["capacity_onshoreA"]
    WindCapB = typeof(WindData["capacity_onshoreB"]) == Float64 ? [WindData["capacity_onshoreB"]] : WindData["capacity_onshoreB"]
    WindCapTotal = WindCapA .+ WindCapB
    replace_nans!(WindCapTotal)

    WindCapOff = typeof(WindData["capacity_offshore"]) == Float64 ? [WindData["capacity_offshore"]] : WindData["capacity_offshore"]
    replace_nans!(WindCapOff)

    PVCapA = typeof(SolarData["capacity_pvplantA"]) == Float64 ? [SolarData["capacity_pvplantA"]] : SolarData["capacity_pvplantA"]
    PVCapB = typeof(SolarData["capacity_pvplantB"]) == Float64 ? [SolarData["capacity_pvplantB"]] : SolarData["capacity_pvplantB"]
    CSPCapA = typeof(SolarData["capacity_cspplantA"]) == Float64 ? [SolarData["capacity_cspplantA"]] : SolarData["capacity_cspplantA"]
    CSPCapB = typeof(SolarData["capacity_cspplantB"]) == Float64 ? [SolarData["capacity_cspplantB"]] : SolarData["capacity_cspplantB"]
    RooftopCap = typeof(SolarData["capacity_pvrooftop"]) == Float64 ? [SolarData["capacity_pvrooftop"]] : SolarData["capacity_pvrooftop"]
    SolarCapTotal = PVCapA .+ PVCapB .+ CSPCapA .+ CSPCapB .+ RooftopCap
    replace_nans!(SolarCapTotal)

    HydroCap = typeof(HydroData["existingcapac"]) == Float64 ? [HydroData["existingcapac"]] : HydroData["existingcapac"]
    replace_nans!(HydroCap)

    # Get the Production data per VRE
    
    # Generating Solar Production Data for each node from the Matlab files
    # Node Number is not save

    PVProductionA = typeof(SolarData["CFtime_pvplantA"]) == Float64 ? [SolarData["CFtime_pvplantA"]] : SolarData["CFtime_pvplantA"]
    PVProductionB = typeof(SolarData["CFtime_pvplantB"]) == Float64 ? [SolarData["CFtime_pvplantB"]] : SolarData["CFtime_pvplantB"]
    CSPProductionA = typeof(SolarData["CFtime_cspplantA"]) == Float64 ? [SolarData["CFtime_cspplantA"]] : SolarData["CFtime_cspplantA"]
    CSPProductionB = typeof(SolarData["CFtime_cspplantB"]) == Float64 ? [SolarData["CFtime_cspplantB"]] : SolarData["CFtime_cspplantB"]
    PVrooftop = typeof(SolarData["CFtime_pvrooftop"]) == Float64 ? [SolarData["CFtime_pvrooftop"]] : SolarData["CFtime_pvrooftop"]
    SolarProductionTotal = PVProductionA .+ PVProductionB .+ CSPProductionA .+ CSPProductionB .+ PVrooftop
    replace_nans!(SolarProductionTotal)

    ## Wind Onshore
    WindProductionA = typeof(WindData["CFtime_windonshoreA"]) == Float64 ? [WindData["CFtime_windonshoreA"]] : WindData["CFtime_windonshoreA"]
    WindProductionB = typeof(WindData["CFtime_windonshoreB"]) == Float64 ? [WindData["CFtime_windonshoreB"]] : WindData["CFtime_windonshoreB"]
    WindProductionOn = WindProductionA .+ WindProductionB
    replace_nans!(WindProductionOn)

    ## Wind Offshore
    WindProductionOff = typeof(WindData["CFtime_windoffshore"]) == Float64 ? [WindData["CFtime_windoffshore"]] : WindData["CFtime_windoffshore"]
    replace_nans!(WindProductionOff)
   
    # demand
    gisdemand = JLD.load(joinpath(inputdata, "SyntheticDemand_$(regionset)_$(sspscenario)-$(sspyear)_$(era_year).jld"), "demand")

    # Transmission distances
    Distances = typeof(TransmissionData["distances"]) == Float64 ? [TransmissionData["distances"]] : TransmissionData["distances"]
    Regionlist = typeof(TransmissionData["regionlist"]) == Float64 ? [TransmissionData["regionlist"]] : TransmissionData["regionlist"]

    
    # Hydro
    n = length(Regionlist)

    existingcapac = permutedims(HydroData["existingcapac"]) .* 1000
    #Monthly inflow
    existinginflow = permutedims(HydroData["existinginflowcf"])
    # Replace NaNs from monthly data with 0's
    replace_nans!(existinginflow)
    # Preallocate array for intermediate result
    avail_inflow = zeros(Float64, 8760, n)           ## Make it variable

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
    reservoirp = [0 0.953 1 0 0.433]                                    ## where does this come from?

    #Hydro capacities
    hydroCap = existingcapac .* reservoirp
    hydroRoRCap = existingcapac - hydroCap
    hydrocapacity = [hydroCap; hydroRoRCap]
    
    #Reservoir and RoR inflow
    hyd_in = avail_inflow .* hydroCap
    hydRoR_in = avail_inflow .* hydroRoRCap

end


function create_node_data()
   
    # number of nodes and technologies 
    n = length(Regionlist)
    t = 9                   # eventually try to find a way to not have it hard coded
    
    # Create a Vector for every technology spanning all nodes (at 5 nodes: 43.800 columns)
    Demand = []
    WindProdOff = []
    WindProdOn = []
    SolarProd = []
    hydro_in = []
    hydroRoR_in = []
    
    for i in 1:n
        Demand$i = sum(gisdemand[:,i], dims = 2)
        append!(Demand, Demand$i)
    end
    for i in 1:n
        WindProdOff$i = sum(WindProductionOff[:,:,i], dims = 2)*1000
        append!(WindProdOff, WindProdOff$i)
    end
    for i in 1:n
        WindProdOn$i = sum(WindProductionOn[:,:,i], dims = 2)*1000
        append!(WindProdOn, WindProdOn$i)
    end
    for i in 1:n
        SolarProd$i = sum(SolarProductionTotal[:,:,i], dims = 2)*1000
        append!(SolarProd, SolarProd$i)
    end
    for i in 1:n
        hyd_in$i = sum(hyd_in[:,i], dims = 2)*1000
        append!(hydro_in, hyd_in$i)
    end
    for i in 1:n
        hydRoR_in$i = sum(hydRoR_in[:,i], dims = 2)*1000
        append!(hydroRoR_in, hydRoR_in$i)
    end
    
    
    All = hcat(Demand, SolarProd, WindProdOn, WindProdOff, hydro_in, hydroRoR_in) # hydro_in, hydroRoR_in
    
    # Generate CSV files for every single node (i.e. a new file every 8760 hours) 
    j = 8760
    x = 1
    
    for m = 1:n
        NodeData = All[x:j*m,:]
        x = x+8760
        Header = ["Demand" "Avail_Sol" "Avail_Wind_On" "Avail_Wind_Off" "Hyd_In" "HydRoR_In"] #"Hyd_In" "HydRoR_In"
        N_m = vcat(Header, NodeData)
        writedlm("N_$m.csv", N_m, ',')
    end
    
    
    # create file for node_specs
    NodeNr = [1:n;]
    Max_Demand = repeat(1:1; outer=[length(Regionlist)])
    nodes_specs_body = hcat(NodeNr, Regionlist, Max_Demand)
    HeaderN = ["Node" "Name" "Max_Demand" ]
    nodes_specs = vcat(HeaderN, nodes_specs_body)
    writedlm("nodes_specs.csv", nodes_specs, ',')
end
    

function gen_capacity()
    # capacity of the remaining non-renewable energy sources fixed at 1000 GW (biomass, nuclear, coal, gas_cc, gas_oc)
    Capacity_else = 1000*1000
    
    # number of nodes and technologies 
    n = length(Regionlist)
    t = 9                   # eventually try to find a way to not have it hard coded
    
    Capacity = zeros(5,1)
    gcap_max = []
            
    # create a vector with the capacity of each energy source per node
    for i in 1:n
        Capacity$i = [sum(WindCapTotal[i,:]), sum(WindCapOff[i,:]), sum(SolarCapTotal[i,:]),  Capacity_else, Capacity_else, Capacity_else, Capacity_else, Capacity_else, sum(hydrocapacity[:,i])].*1000
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
    Header = ["node" "gen_tech" "gcap_min" "gcap_max(GW)"]
    generation_capacity = vcat(Header, Gen_capac)
    
    writedlm("gen_capacity.csv", generation_capacity, ',')
    
end


function sto_capacity()
    n = length(Regionlist)
    st = 1          
    Max_capacity = 1000000000
        
    s = repeat([1:st;], n)
    node = repeat(1:n, inner=st)
    scap_min = zeros(n)
    scap_max = fill(Max_capacity, (n*st, 1))
    sto_cap = hcat(s, node, scap_min, scap_max)
    Header = ["s" "node" "scap_min" "scap_max"]
    sto_capacity = vcat(Header, sto_cap)
    
    writedlm("sto_capacity.csv", sto_capacity, ',')
end


function transmission()
    Dist = Distances[tril!(trues(size(Distances)), -1)]

    n = length(Regionlist)
    
    # get the combinations of the pair of nodes
    Nodes = [1:n;]
    Pairs = collect(combinations(Nodes, 2))
    L_ind = [1:length(Pairs);]
    
    # Create Matrix for remaining parameters of transmission
    
    TransData = [
       # :cost       :converter_cost         :M      :C      :B      :efficiency         :lifetime       :tcap_min       :tcap_max
        400          150000                  0.02    0       1       0.95                40              0               1000000000
    ]
    
    # Fill the matrix according to the number of pairs in the system
    TransDatafull = repeat(TransData; outer=[length(Pairs)])
    
    # combine the tables, add header and create CSV file
    TransmissionData = hcat(L_ind, Pairs, Dist, TransDatafull)
    Header = ["L_ind" "line" "dist"  "cost" "converter_cost" "M" "C"  "B" "efficiency" "lifetime" "tcap_min" "tcap_max"]
    TransmissionTable = vcat(Header, TransmissionData)
    writedlm("transmission.csv", TransmissionTable, ',')
end


function generate_tech_tables()

    techtable = [
        :name          :g       :investment_cost    :fixedOM    :varOM      :fuel_cost      :efficiency     :emissions      :lifetime       :r_minus    :r_plus
        :wind_on        1       1127                35          0           0               1               0               25              1           1
        :wind_off       2       2290                80          0           0               1               0               25              1           1
        :solar          3       480                 26          0           0               1               0               25              1           1
        :biomass        4       2076                100         0           7               0.448           0.39            30              1           1
        :nuclear        5       5000                150         3           3               0.34            0               60              0.05        0.05
        :coal           6       1300                25          6           8               0.466           0.34            40              0.15        0.15
        :gas_cc         7       800                 20          4           22              0.615           0.20            30              0.3         0.3
        :gas_oc         8       400                 15          3           22              0.395           0.20            30              1           1
    ]
    
    writedlm("gen_technology.csv", techtable, ',')
        
    # Hydro
    hydro_tech = [
        :name          :hydro_tech       :investment_cost    :fixedOM    :varOM      :fuel_cost      :efficiency     :emissions      :lifetime       :r_minus    :r_plus
        :hydropower     1                 4000                200          0           0               1               0               60              1           1
    ]
    
    writedlm("hydro_technology.csv", hydro_tech, ',')
        
    # Storage
    storage = [
        :s          :xi      :c    :cost    :lifetime        
        1           0.86     15     185     20
    ]
    writedlm("storage.csv", storage, ',')
   
end


function create_data_sets(inputdata, regionset, sspscenario, sspyear, era_year)
    get_data(inputdata, regionset, sspscenario, sspyear, era_year)
    create_node_data()
    gen_capacity()
    sto_capacity()
    transmission()
    generate_tech_tables()
end