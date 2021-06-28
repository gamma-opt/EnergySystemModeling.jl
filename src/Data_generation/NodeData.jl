
using MAT, DelimitedFiles, JLD

# path to your set output folder of GlobalEnergyGIS data
inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

# Name of the regionset you defined in Inputdata.jl
regionset = "EuropeTest"

SolarData =matread(joinpath(inputdata, "GISdata_solar2018_$regionset.mat"))

WindData = matread(joinpath(inputdata, "GISdata_wind2018_$regionset.mat"))

HydroData = matread(joinpath(inputdata, "GISdata_hydro_$regionset.mat"))


# Generating Solar Production Data for each node from the Matlab files
# Node Number is not save

SolarProductionA = typeof(SolarData["CFtime_pvplantA"]) == Float64 ? [SolarData["CFtime_pvplantA"]] : SolarData["CFtime_pvplantA"]
SolarProductionB = typeof(SolarData["CFtime_pvplantB"]) == Float64 ? [SolarData["CFtime_pvplantB"]] : SolarData["CFtime_pvplantB"]
SolarProductionTotal = SolarProductionA .+ SolarProductionB
SolarProductionTotal[isnan.(SolarProductionTotal)].= 0


## Wind Onshore
WindProductionA = typeof(WindData["CFtime_windonshoreA"]) == Float64 ? [WindData["CFtime_windonshoreA"]] : WindData["CFtime_windonshoreA"]
WindProductionB = typeof(WindData["CFtime_windonshoreB"]) == Float64 ? [WindData["CFtime_windonshoreB"]] : WindData["CFtime_windonshoreB"]
WindProductionOn = WindProductionA .+ WindProductionB
WindProductionOn[isnan.(WindProductionOn)].= 0


## Wind Offshore

WindProductionOff = typeof(WindData["CFtime_windoffshore"]) == Float64 ? [WindData["CFtime_windoffshore"]] : WindData["CFtime_windoffshore"]
WindProductionOff[isnan.(WindProductionOff)].= 0


## Hydro
WindProductionOff = typeof(WindData["CFtime_windoffshore"]) == Float64 ? [WindData["CFtime_windoffshore"]] : WindData["CFtime_windoffshore"]
WindProductionOff[isnan.(WindProductionOff)].= 0


## demand
gisdemand = JLD.load(joinpath(inputdata, "SyntheticDemand_$(regionset)_ssp2-26-2050_2018.jld"), "demand")

# Get the different regions created in GlobalEnergyGIS
TransmissionData = matread(joinpath(inputdata, "distances_$regionset.mat"))
Regionlist = typeof(TransmissionData["regionlist"]) == Float64 ? [TransmissionData["regionlist"]] : TransmissionData["regionlist"]

# number of nodes and technologies 
n = length(Regionlist)
t = 8

# Create a Vector for every technology spanning all nodes (at 5 nodes: 43.800 columns)
Demand = []
WindProdOff = []
WindProdOn = []
SolarProd = []

for i in 1:n
    Demand$i = gisdemand[:,i]
    append!(Demand, Demand$i)
end

for i in 1:n
    WindProdOff$i = WindProductionOff[:,1,i]
    append!(WindProdOff, WindProdOff$i)
end

for i in 1:n
    WindProdOn$i = WindProductionOn[:,1,i]
    append!(WindProdOn, WindProdOn$i)
end

for i in 1:n
    SolarProd$i = SolarProductionTotal[:,1,i]
    append!(SolarProd, SolarProd$i)
end

All = hcat(Demand, SolarProd, WindProdOn, WindProdOff)

# Generate CSV files for every single node (i.e. a new file every 8760 hours) 
j = 8760
x = 1

for m = 1:n
    NodeData = All[x:j*m,:]
    x = x+8760
    Header = ["Demand" "Solar" "WindOn" "WindOff"]
    N_$m = vcat(Header, NodeData)
    writedlm("N_$m.csv", NodeData, ',')
end


# create file for node_specs
NodeNr = [1:n;]
Max_Demand = repeat(1:1; outer=[length(Regionlist)])
nodes_specs_body = hcat(NodeNr, Regionlist, Max_Demand)
HeaderN = ["Node" "Name" "Max_Demand" ]
nodes_specs = vcat(HeaderN, nodes_specs_body)
writedlm("nodes_specs.csv", nodes_specs, ',')