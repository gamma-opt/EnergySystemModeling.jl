using DataFrames: Matrix
using MAT, DelimitedFiles, JLD, JSON, CSV, LinearAlgebra, Combinatorics, GlobalEnergyGIS, DataFrames

include(joinpath("src", "data_generation", "data_generation.jl"))


regionset = "Regions"       
sspscenario_input = "ssp2-26"
sspyear_input = 2050
era_year_input = 2018

# Path to the output folder of GlobalEnergyGIS
inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

# Write the countries for which you want to create the dataset. 
Regions = ["Germany", "Austria", "Switzerland"]

# use either Country_names or Regions
Dataset =  get_countries(Regions)

create_data_sets(inputdata, regionset, sspscenario_input, sspyear_input, era_year_input, Dataset)



