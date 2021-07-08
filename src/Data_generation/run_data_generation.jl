using MAT, DelimitedFiles, JLD, JSON, CSV, LinearAlgebra, Combinatorics, GlobalEnergyGIS, DataFrames


regionset = "EU27"       
sspscenario_input = "ssp2-26"
sspyear_input = 2050
era_year_input = 2018

# Write the countries for which you want to create the dataset
Country_names = ["Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia",
                "Finland", "France", "Germany", "Hungary", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg",
                "Malta", "Netherlands", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden"]

# Name of your Dataset



# Path to the output folder of GlobalEnergyGIS
inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"

Dataset =  receive_countries(Country_names) 

create_data_sets(inputdata, regionset, sspscenario_input, sspyear_input, era_year_input, Dataset)



