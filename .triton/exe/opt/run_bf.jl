using Logging
using EnergySystemModeling
using JLD2
using JuMP
using Gurobi
using Plots
using StatsPlots
using Dates
using DelimitedFiles

@info string("Starting instance ",ARGS[1]," @ ",now()) 

WRKDIR = "/scratch/work/condeil1/EnergySystemModeling.jl/examples"
cd(WRKDIR)

@info "Loading parameters"
constants_dir = "constants"
constants_path = joinpath(WRKDIR,constants_dir)
structure_dir = "8nodes"
structures_path = joinpath(WRKDIR,"structures",structure_dir)
instances_path = joinpath(structures_path,"instances/.big_files")
instance_clust = ARGS[1]
instance_path_clust = joinpath(instances_path,instance_clust)
instance_ftr = "08n8760h_ftr"
instance_path_ftr = joinpath(instances_path,instance_ftr)

@info "Creating output directories"
output_path = joinpath(instance_path_clust,"output_tri_bf75")
output_ftr_path = joinpath(instance_path_ftr,"output_tri")

results_dir = "results"
plots_dir = "plots"
plots_output_path = joinpath(output_path,plots_dir)
csv_dir = "csv"
model_dir = joinpath(".big_files","model")
params_path_ftr = joinpath(instances_path,instance_ftr,".big_files","parameters")

mkpath(output_path)
mkpath(joinpath(output_path,results_dir))

## Other possible folders
# mkpath(joinpath(output_path,plots_dir))
# mkpath(joinpath(output_path,plots_dir,"pdf"))
# mkpath(joinpath(output_path,plots_dir,"png"))
# mkpath(joinpath(output_path,csv_dir))
# mkpath(joinpath(output_path,model_dir))

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

@info "Extracting clust. variables (FTR)"
variables = load(joinpath(output_ftr_path, results_dir, "variables.jld2"));

@info "Creating the energy system model"
(model, VariablesDict, ObjectivesDict) = EnergySystemModel(parameters, specs)


@info "Updating lower bounds (brownfield)"
# Updating lower bounds (BROWNFIELD)
cap_percentual = 1

for g in parameters.G, n in parameters.N
    set_lower_bound(model[:p̄_gn][g,n], cap_percentual.*variables["p̄_gn"][g,n])
end
# if specs.transmission
#     for l in parameters.L
#         set_lower_bound(model[:f̄_l][l], cap_percentual*variables["f̄_l"][l])
#     end
# end
# if specs.storage
#     for s in parameters.S, n in parameters.N
#         set_lower_bound(model[:b̄_sn][s,n], cap_percentual.*variables["b̄_sn"][s,n])
#     end
# end
# if specs.hydro
#     for h in parameters.H, n in parameters.N
#         set_lower_bound(model[:h̄_hn][h,n], cap_percentual.*variables["h̄_hn"][h,n])
#     end
# end

@info "Optimizing the model"
optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60*60*72,
                                      "LogFile" => joinpath(output_path, "gurobi.log"))
set_optimizer(model, optimizer)
set_optimizer_attributes(model, "NumericFocus" => 1)
set_optimizer_attributes(model, "ScaleFlag" => 2)

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

## Write a out file
# open(joinpath(instance_path_clust,"out1.txt"), "w") do io
#         writedlm(io, "")
# end

@info string("Ended instance ",ARGS[1]," @ ",now())