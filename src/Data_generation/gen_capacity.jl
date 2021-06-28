using MAT, DelimitedFiles, JLD

    # Path to the output folder for the data of GlobalEnergyGIS
    inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

    SolarData = matread(joinpath(inputdata, "GISdata_solar2018_EuropeSmall.mat"))

    WindData = matread(joinpath(inputdata, "GISdata_wind2018_EuropeSmall.mat"))

    HydroData = matread(joinpath(inputdata, "GISdata_hydro_EuropeSmall.mat"))

    # Get the different capacities per VRE source from GlobalEnergyGIS files
    WindCapA = typeof(WindData["capacity_onshoreA"]) == Float64 ? [WindData["capacity_onshoreA"]] : WindData["capacity_onshoreA"]
    WindCapB = typeof(WindData["capacity_onshoreB"]) == Float64 ? [WindData["capacity_onshoreB"]] : WindData["capacity_onshoreB"]
    WindCapTotal = WindCapA .+ WindCapB
    WindCapTotal[isnan.(WindCapTotal)].= 0

    WindCapOff = typeof(WindData["capacity_offshore"]) == Float64 ? [WindData["capacity_offshore"]] : WindData["capacity_offshore"]
    WindCapOff[isnan.(WindCapTotal)].= 0

    SolarCapA = typeof(SolarData["capacity_pvplantA"]) == Float64 ? [SolarData["capacity_pvplantA"]] : SolarData["capacity_pvplantA"]
    SolarCapB = typeof(SolarData["capacity_pvplantB"]) == Float64 ? [SolarData["capacity_pvplantB"]] : SolarData["capacity_pvplantB"]
    SolarCapTotal = SolarCapA .+ SolarCapB
    SolarCapTotal[isnan.(SolarCapTotal)].= 0

    HydroCap = typeof(HydroData["existingcapac"]) == Float64 ? [HydroData["existingcapac"]] : HydroData["existingcapac"]
    HydroCap[isnan.(HydroCap)].= 0

    # capacity of the remaining non-renewable energy sources fixed at 1000 GW
    Capacity_else = 1000

    # number of nodes and technologies 
    n = 5
    t = 8

    Capacity = zeros(5,1)
    gcap_max = []
    
    # create a vector with the capacity of each energy source per node
    for i in 1:n
        Capacity$i = [WindCapTotal[i,1], WindCapOff[i,1], SolarCapTotal[i,1], HydroCap[i], Capacity_else, Capacity_else, Capacity_else, Capacity_else]
        append!(gcap_max, Capacity$i)
    end
    
    gen_capacity = zeros(n*t, 4)
    
    node = repeat(1:n, inner=t)
    gen_tech = repeat([1:t;], n)
    gcap_min = fill(0, (40,1))

  
    Gen_capac = hcat(node, gen_tech, gcap_min, gcap_max)
    Header = ["node" "gen_tech" "gcap_min" "gcap_max(GW)"]
    gen_capacity = vcat(Header, Gen_capac)

    writedlm("gen_capacity2.csv", gen_capacity, ',')

