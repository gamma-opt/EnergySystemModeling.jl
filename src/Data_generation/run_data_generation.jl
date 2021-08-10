using Combinatorics: length
using DataFrames: Matrix
using MAT, DelimitedFiles, JLD, JSON, CSV, LinearAlgebra, Combinatorics, GlobalEnergyGIS, DataFrames, Statistics

include("data_generation.jl")

# name the instance that you're creating
sspscenario_input = "ssp2-26"
sspyear_input = 2050
era_year_input = 2018

# select the folder in which you want to write the files, "folder" and "subfolder" should not change
folder = "examples"
subfolder = "structures"
instance = "5_regions"

# Path to the output folder of GlobalEnergyGIS
inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

# Write the countries/regions for which you want to create the dataset. Use either country Names or predifined regions: Nordics, Eastern, Western, Central, Mediterranian
Regions = ["Nordics", "Eastern", "Western", "Central", "Mediterranian"]

# Desired time period in hours
T = 8760
# number of technologies
t = 8       # 8 since hydro is not included here 

Dataset = get_countries(Regions)

create_data_sets(inputdata, sspscenario_input, sspyear_input, era_year_input, Dataset, folder, subfolder, instance, T, t)



