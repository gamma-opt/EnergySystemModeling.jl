using Combinatorics: length
using DataFrames: Matrix
using MAT, DelimitedFiles, JLD, JSON, CSV, LinearAlgebra, Combinatorics, GlobalEnergyGIS, DataFrames, Statistics

include("data_generation.jl")

# name the instance that you're creating
regionset = "Regions"       
sspscenario_input = "ssp2-26"
sspyear_input = 2050
era_year_input = 2018

# select the folder in which you want to write the files, "folder" and "subfolder" should not change
folder = "examples"
subfolder = "structures"
instance_name = "5_regions"

# Path to the output folder of GlobalEnergyGIS
inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

# Write the countries/regions for which you want to create the dataset. Use either country Names or predifined regions: Nordics, Eastern, Western, Central, Mediterranian
Regions = ["Nordics", "Mediterranian", "Western", "Eastern", "Central"]

Dataset = get_countries(Regions)

create_data_sets(inputdata, regionset, sspscenario_input, sspyear_input, era_year_input, Dataset)



