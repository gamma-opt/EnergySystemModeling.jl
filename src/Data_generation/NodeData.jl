
using MAT, DelimitedFiles, JLD

# path to your set output folder of GlobalEnergyGIS data
inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

SolarData =matread(joinpath(inputdata, "GISdata_solar2018_EuropeSmall.mat"))

WindData = matread(joinpath(inputdata, "GISdata_wind2018_EuropeSmall.mat"))

HydroData = matread(joinpath(inputdata, "GISdata_hydro_EuropeSmall.mat"))


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
gisdemand = JLD.load(joinpath(inputdata, "SyntheticDemand_EuropeSmall_ssp2-26-2050_2018.jld"), "demand")


n = 5
t = 8

Demand = []
WindProdOff = []
WindProdOn = []
SolarProd = []

# create a vector with the capacity of each energy source per node
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


j = 8760
x = 1

for m = 1:n
    NodeData = All[x:j*m,:]
    x = x+8760
    Header = ["Demand" "Solar" "WindOn" "WindOff"]
    N_$m = vcat(Header, NodeData)
    writedlm("N_$m.csv", NodeData, ',')
end


