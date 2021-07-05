using MAT, DelimitedFiles, JLD, JSON, CSV

regionset = "EuropeTest"
sspscenario = "ssp2-26"
sspyear = 2050
era_year = 2018

# Path to the output folder of GlobalEnergyGIS
inputdata = "D:\\Eigene Dateien\\Studium\\Master\\RA\\Copernicus\\output"


create_data_sets(inputdata, regionset, sspscenario, sspyear, era_year)
