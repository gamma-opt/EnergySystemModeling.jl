# Data generation

## Run_data_generation.jl

Users set the framework for their desired regionsets and generate the regionset.
The instance is then created using the previously defined framework and dataset by: 
```@example create_data_sets
create_data_sets(inputdata, regionset, sspscenario_input, sspyear_input, era_year_input, Dataset, folder, subfolder, instance_name)
```

## Data_generation.jl

### Get country list
```@example get_countries
get_countries(Regions)
```
Takes the countries or regions specified by the user and turns them into a structure which GlobalEnergyGIS can recognize and on base of which it creates datasets of renewable energies.

```@example get_countries
get_countries()
```
If the user does not specify any set of countries or regions, he will get an instance for whole Europe.

### Generate datasets and CSV files
```@example create_data_sets
create_data_sets(inputdata, sspscenario_input, sspyear_input, era_year_input, Dataset, folder, subfolder, instance, T, t)
```

Uses GlobalEnergyGIS to generate renewable energy input data based on the sspscenario, sspyear and base year defined by the user. It then reads the generated data, converts it and generates different CSV files with data that can be used by EnergySystemModeling.jl.