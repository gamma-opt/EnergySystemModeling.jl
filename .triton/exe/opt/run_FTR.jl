using Logging
using EnergySystemModeling
using JLD2
using JuMP
using Gurobi
using Plots
using StatsPlots
using Dates
using DelimitedFiles

@info string("Starting FTR instance @ ",now())

WRKDIR = "/scratch/work/condeil1/EnergySystemModeling.jl/examples"
cd(WRKDIR)

@info "Loading parameters"
constants_dir = "constants"
constants_path = joinpath(WRKDIR,constants_dir)
structure_dir = "8nodes"
structures_path = joinpath(WRKDIR,"structures",structure_dir)
instances_path = joinpath(structures_path,"instances/.big_files")
instance_clust = "08n8760h_ftr"
instance_path_clust = joinpath(instances_path,instance_clust)
instance_ftr = "08n8760h_ftr"

@info "Creating output directories"
output_path = joinpath(instance_path_clust,"output_tri")
results_dir = "results"
plots_dir = "plots"
plots_output_path = joinpath(output_path,plots_dir)
csv_dir = "csv"
model_dir = joinpath(".big_files","model")
params_path_ftr = joinpath(instances_path,instance_ftr,".big_files","parameters")

mkpath(output_path)
mkpath(joinpath(output_path,results_dir))
mkpath(joinpath(output_path,plots_dir))
mkpath(joinpath(output_path,plots_dir,"pdf"))
mkpath(joinpath(output_path,plots_dir,"png"))
mkpath(joinpath(output_path,csv_dir))
mkpath(joinpath(output_path,model_dir))

@info "Creating specifications"
# Updating the parameters to the clustering instance
parameters = change_time_parameters(params_path_ftr, instance_path_clust)

# Model specifications
specs = Specs(
        transmission=true,
        renewable_target=true,
        carbon_cap=true,
        nuclear_limit=false,
        storage=true,
        ramping=true,
        voltage_angles=false,
        hydro=true,
        hydro_simple=false
)

## Plotting specifications
Plots_specs = Dict{String,Bool}()
# Plotting part 1: Objective function values
Plots_specs["p1"] = false
# Plotting part 2: Dispatch and storage
Plots_specs["p2"] = false
# Plotting part 3: Storage capacity
Plots_specs["p3"] = false
# Plotting part 4: Generation technologies dispatch levels
Plots_specs["p4"] = false
# Plotting part 5: Generation capacities (including hydro)
Plots_specs["p5"] = false
# Plotting part 6: Consolidated dispatch vs demand
Plots_specs["p6"] = false
# Plotting part 7: Transmission flow (per line)
Plots_specs["p7"] = false
# Plotting part 8: Transmission capacities
Plots_specs["p8"] = false
# Plotting part 9: Consolidated transmission flow
Plots_specs["p9"] = false
# Plotting part 10: Lost of load
Plots_specs["p10"] = false
# Plotting in development
Plots_specs["p11"] = false

@info "Creating the energy system model"
(model, VariablesDict, ObjectivesDict) = EnergySystemModel(parameters, specs)
# write_to_file(model,joinpath(output_path, model_dir,"model.mps"))

@info "Optimizing the model"
optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60*60*72,
                                      "LogFile" => joinpath(output_path, "gurobi.log"))
set_optimizer(model, optimizer)
set_optimizer_attributes(model, "NumericFocus" => 1)
set_optimizer_attributes(model, "ScaleFlag" => 2)

# Scaling: ("ScaleFlag" => [-1,3])
# Numeric focus: ("NumericFocus" => [0,3])
# set_optimizer_attributes(model, "Method" => 0)
# set_optimizer_attributes(model, "Crossover" => 1)

optimize!(model)

@info "Extracting results"
variables = JuMPVar(model, VariablesDict)
objectives = JuMPObj(model, ObjectivesDict)
expressions = Expressions(parameters, specs, variables)

# @info "Saving results (JSON)"
# save_json(specs, joinpath(output_path, results_dir, "specs.json"))
# save_json(variables, joinpath(output_path, results_dir, "variables.json"))
# save_json(objectives, joinpath(output_path, results_dir, "objectives.json"))
# save_json(expressions, joinpath(output_path, results_dir, "expressions.json"))

@info "Saving results (JLD2)"
JLD2.save(joinpath(output_path, results_dir, "variables.jld2"), variables; compress = true)
JLD2.save(joinpath(output_path, results_dir, "objectives.jld2"), objectives; compress = true)
JLD2.save(joinpath(output_path, results_dir, "expressions.jld2"), expressions; compress = true)

# @info "Perform plotting"
# perform_plotting(Plots_specs, parameters, variables, objectives, expressions, plots_output_path)

## Write a out file
# open(joinpath(instance_path_clust,"out1.txt"), "w") do io
#         writedlm(io, "")
# end

@info string("Ended FTR instance @ ",now())