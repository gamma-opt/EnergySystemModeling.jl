using MAT, DelimitedFiles, JLD, JSON, CSV, LinearAlgebra, Combinatorics, GlobalEnergyGIS

regionset = "EU27"       # If changed here, you need to change it in the saveregions function of data_generation.jl
sspscenario = "ssp2-26"
sspyear = 2050
era_year = 2018

Dataset = "Dataset_EU27"   # If changed here, you need to change it in the saveregions function of data_generation.jl

# Path to the output folder of GlobalEnergyGIS
inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"


create_data_sets(inputdata, regionset, sspscenario, sspyear, era_year, Dataset)

## List of countries with their SI code: https://gist.github.com/tadast/8827699#file-countries_codes_and_coordinates-csv


