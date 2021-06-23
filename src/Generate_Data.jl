
using GlobalEnergyGIS, MAT, DelimitedFiles, JLD

create_scenario_datasets("SSP2", 2050)

EuropeSmall = [
    "Nordics"           GADM("Norway", "Sweden", "Finland", "Denmark")
    "Mediterranian"     GADM("Spain", "Portugal", "Italy", "Croatia", "Greece", "Bosnia and Herzegovina", "Slowenia", "Serbia", "Kosovo", "Albania", "Macedonia")
    "Western"           GADM("France", "Belgium", "Netherland", "Luxemburg")
    "Central"           GADM("Germany", "Austria", "Switzerland")    
    "Eastern"           GADM("Poland", "Czechia", "Slovakia", "Hungary")
]

saveregions("EuropeSmall", EuropeSmall)

makedistances("EuropeSmall")

createmaps("EuropeSmall")




inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

SolarData =matread(joinpath(inputdata, "GISdata_solar2018_EuropeTest.mat"))

WindData = matread(joinpath(inputdata, "GISdata_wind2018_EuropeTest.mat"))

HydroData = matread(joinpath(inputdata, "GISdata_hydro_EuropeTest.mat"))


# Generating Solar Production Data for each node from the Matlab files
# Node Number is not save


SolarProductionA = typeof(SolarData["CFtime_pvplantA"]) == Float64 ? [SolarData["CFtime_pvplantA"]] : SolarData["CFtime_pvplantA"]
SolarProductionB = typeof(SolarData["CFtime_pvplantB"]) == Float64 ? [SolarData["CFtime_pvplantB"]] : SolarData["CFtime_pvplantB"]
SolarProductionTotal = SolarProductionA .+ SolarProductionB
SolarProductionTotal[isnan.(SolarProductionTotal)].= 0

#First row of every "block" / matrix

# Nordics (node 1)
SolarProductionNordics = SolarProductionTotal[:,1,1]
writedlm("SolarProductionNordics.csv", SolarProductionNordics, ',')

# Mediterranian (node 2)
SolarProductionMed = SolarProductionTotal[:,1,2]
writedlm("SolarProductionMed.csv", SolarProductionMed, ',')

# Western (node 3)
SolarProductionWest = SolarProductionTotal[:,1,3]
writedlm("SolarProductionWest.csv", SolarProductionWest, ',')

# Central (node 4)
SolarProductionCen = SolarProductionTotal[:,1,4]
writedlm("SolarProductionCen.csv", SolarProductionCen, ',')

# Eastern (node 5)
SolarProductionEas = SolarProductionTotal[:,1,5]
writedlm("SolarProductionEas.csv", SolarProductionEas, ',')



## Wind Onshore
WindProductionA = typeof(WindData["CFtime_windonshoreA"]) == Float64 ? [WindData["CFtime_windonshoreA"]] : WindData["CFtime_windonshoreA"]
WindProductionB = typeof(WindData["CFtime_windonshoreB"]) == Float64 ? [WindData["CFtime_windonshoreB"]] : WindData["CFtime_windonshoreB"]
WindProductionTotal = WindProductionA .+ WindProductionB
WindProductionTotal[isnan.(WindProductionTotal)].= 0

# Nordics (node 1)
WindProductionNordics = WindProductionTotal[:,1,1]
writedlm("WindProductionNordics.csv", WindProductionNordics, ',')

# Mediterranian (node 2)
WindProductionMed = WindProductionTotal[:,1,2]
writedlm("WindProductionMed.csv", WindProductionMed, ',')

# Western (node 3)
WindProductionWest = WindProductionTotal[:,1,3]
writedlm("WindProductionWest.csv", WindProductionWest, ',')

# Central (node 4)
WindProductionCen = WindProductionTotal[:,1,4]
writedlm("WindProductionCen.csv", WindProductionCen, ',')

# Eastern (node 5)
WindProductionEas = WindProductionTotal[:,1,5]
writedlm("WindProductionEas.csv", WindProductionEas, ',')



## Wind Offshore

WindProductionOff = typeof(WindData["CFtime_windoffshore"]) == Float64 ? [WindData["CFtime_windoffshore"]] : WindData["CFtime_windoffshore"]
WindProductionOff[isnan.(WindProductionOff)].= 0


# Nordics (node 1)
WindProductionNordicsOff = WindProductionOff[:,1,1]
writedlm("WindProductionNordicsOff.csv", WindProductionNordicsOff, ',')

# Mediterranian (node 2)
WindProductionMedOff = WindProductionOff[:,1,2]
writedlm("WindProductionMedOff.csv", WindProductionMedOff, ',')

# Western (node 3)
WindProductionWestOff = WindProductionOff[:,1,3]
writedlm("WindProductionWestOff.csv", WindProductionWestOff, ',')

# Central (node 4)
WindProductionCenOff = WindProductionOff[:,1,4]
writedlm("WindProductionCenOff.csv", WindProductionCenOff, ',')

# Eastern (node 5)
WindProductionEasOff = WindProductionOff[:,1,5]
writedlm("WindProductionEasOff.csv", WindProductionEasOff, ',')


## Hydro
WindProductionOff = typeof(WindData["CFtime_windoffshore"]) == Float64 ? [WindData["CFtime_windoffshore"]] : WindData["CFtime_windoffshore"]
WindProductionOff[isnan.(WindProductionOff)].= 0


## demand

gisdemand = JLD.load(joinpath(inputdata, "SyntheticDemand_EuropeTest_ssp2-34-2050_2018.jld"), "demand")
DemandNordics = gisdemand[:,1]

DemandMed = gisdemand[:,2]

DemandWest = gisdemand[:,3]

DemandCen = gisdemand[:,4]

DemandEas = gisdemand[:,5]


# Concatenate the data for each node

N1_0 = hcat(DemandNordics, SolarProductionNordics, WindProductionNordics, WindProductionNordicsOff)
Header = ["Demand" "Solar" "WindOn" "WindOff"]
N1 = vcat(Header, N1_0)
writedlm("N1.csv", N1, ',')

N2_0 = hcat(DemandMed, SolarProductionMed, WindProductionMed, WindProductionMedOff)
Header = ["Demand" "Solar" "WindOn" "WindOff"]
N2 = vcat(Header, N2_0)
writedlm("N2.csv", N2, ',')

N3_0 = hcat(DemandWest, SolarProductionWest, WindProductionWest, WindProductionWestOff)
Header = ["Demand" "Solar" "WindOn" "WindOff"]
N3 = vcat(Header, N3_0)
writedlm("N3.csv", N3, ',')

N4_0 = hcat(DemandCen, SolarProductionCen, WindProductionCen, WindProductionCenOff)
Header = ["Demand" "Solar" "WindOn" "WindOff"]
N4 = vcat(Header, N4_0)
writedlm("N4.csv", N4, ',')

N5_0 = hcat(DemandEas, SolarProductionEas, WindProductionEas, WindProductionEasOff)
Header = ["Demand" "Solar" "WindOn" "WindOff"]
N5 = vcat(Header, N5_0)
writedlm("N5.csv", N5, ',')
