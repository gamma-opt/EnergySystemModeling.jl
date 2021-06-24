using MAT, DelimitedFiles, JLD


inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

SolarData =matread(joinpath(inputdata, "GISdata_solar2018_EuropeSmall.mat"))

WindData = matread(joinpath(inputdata, "GISdata_wind2018_EuropeSmall.mat"))

HydroData = matread(joinpath(inputdata, "GISdata_hydro_EuropeSmall.mat"))



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

Capacity_else = 1000000000
Gen_capacity = [WindCapTotal[1,1], WindCapOff[1,1], SolarCapTotal[1,1], HydroCap[1], ]
