function replace_nans!(array::Array{Float64, N}) where N
    for i in eachindex(array)
        if isnan(array[i])
            array[i] = zero(1)
        end
    end
end

T = 8760
n = length(Regionlist)
t = 9   

function get_data(inputdata, regionset, sspscenario, sspyear, era_year)

    # Name of the regionset you defined in Inputdata.jl
    SolarData = matread(joinpath(inputdata, "GISdata_solar$(era_year)_$regionset.mat"))
    WindData = matread(joinpath(inputdata, "GISdata_wind$(era_year)_$regionset.mat"))
    HydroData = matread(joinpath(inputdata, "GISdata_hydro_$regionset.mat"))
    TransmissionData = matread(joinpath(inputdata, "distances_$regionset.mat"))

    
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
    replace_nans!(WindCapTotal)

    WindCapOff = typeof(WindData["capacity_offshore"]) == Float64 ? [WindData["capacity_offshore"]] : WindData["capacity_offshore"]
    replace_nans!(WindCapOff)

    ## Wind Onshore
    CFtime_windonshoreA = typeof(WindData["CFtime_windonshoreA"]) == Float64 ? [WindData["CFtime_windonshoreA"]] : WindData["CFtime_windonshoreA"]
    CFtime_windonshoreB = typeof(WindData["CFtime_windonshoreB"]) == Float64 ? [WindData["CFtime_windonshoreB"]] : WindData["CFtime_windonshoreB"]
    CFtime_windonshore = CFtime_windonshoreA .+ CFtime_windonshoreB
    replace_nans!(CFtime_windonshore)

    ## Wind Offshore
    CFtime_windoffshore = typeof(WindData["CFtime_windoffshore"]) == Float64 ? [WindData["CFtime_windoffshore"]] : WindData["CFtime_windoffshore"]
    replace_nans!(CFtime_windoffshore)

      
    #Preallocate arrays for results
    avail_wind_offshore = zeros(Float64, T, n)
    avail_wind_onshoreA = zeros(Float64, T, n)
    avail_wind_onshoreB = zeros(Float64, T, n)

    #Calculate absolute availability values and sum up the different classes for each wind power type
    for i in 1:T
        avail_wind_offshore[i,:,:] = permutedims(sum(WindCapOff .* CFtime_windoffshore[i,:,:], dims=2))
        avail_wind_onshoreA[i,:,:] = permutedims(sum(WindCapA .* CFtime_windonshoreA[i,:,:], dims=2))
        avail_wind_onshoreB[i,:,:] = permutedims(sum(WindCapB .* CFtime_windonshoreB[i,:,:], dims=2))
    end

    avail_wind_on = (avail_wind_onshoreA + avail_wind_onshoreB) ./ capacity_wind_on
    avail_wind_off = avail_wind_offshore ./ capacity_wind_off
    
      
    # demand
    gisdemand = JLD.load(joinpath(inputdata, "SyntheticDemand_$(regionset)_$(sspscenario)-$(sspyear)_$(era_year).jld"), "demand")

    # Transmission distances
    Distances = typeof(TransmissionData["distances"]) == Float64 ? [TransmissionData["distances"]] : TransmissionData["distances"]
    Regionlist = typeof(TransmissionData["regionlist"]) == Float64 ? [TransmissionData["regionlist"]] : TransmissionData["regionlist"]

    # Hydro
    HydroCap = typeof(HydroData["existingcapac"]) == Float64 ? [HydroData["existingcapac"]] : HydroData["existingcapac"]
    replace_nans!(HydroCap)
   
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
    # Create a Vector for every technology spanning all nodes (at 5 nodes: 43.800 columns)
    Demand = []
    WindOff = []
    WindOn = []
    Solar = []
    hydro_in = []
    hydroRoR_in = []
    
    for i in 1:n
        Demand$i = sum(gisdemand[:,i], dims = 2)
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
    
    # Generate CSV files for every single node (i.e. a new file every 8760 hours) 
    
    x = 1
    
    for i = 1:n
        NodeData = All[x:T*i,:]
        x = x+T
        Header = ["Demand" "Avail_Sol" "Avail_Wind_On" "Avail_Wind_Off" "Hyd_In" "HydRoR_In"] #"Hyd_In" "HydRoR_In"
        N_m = vcat(Header, NodeData)
        writedlm("N_$i.csv", N_m, ',')
       # CSV.write(joinpath(DataInput_path, "nodes", "$i.csv"), nodedata)
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