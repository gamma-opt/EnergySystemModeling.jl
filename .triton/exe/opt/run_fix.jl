using Logging
using EnergySystemModeling
using JLD2
using JuMP
using Gurobi
using Dates
using DelimitedFiles

# Define optimisation time limit
opt_tlim = 60*60*72

@info string("Starting instance (fixed capacity) ",ARGS[1]," @ ",now())

# Instance including clustering info
instance_clust = ARGS[1]                # Clustering instance with info
clust_method = ARGS[2]                  # Clustering method ("","dc"/"dc_new"/"dc_nosun"/"day"/etc.)
nosun = parse(Bool,ARGS[3])             # Boolean regarding solar availability

WRKDIR = "/scratch/work/condeil1/EnergySystemModeling.jl/examples"
cd(WRKDIR)

@info "Loading parameters"
constants_dir = "constants"
constants_path = joinpath(WRKDIR,constants_dir)
structure_dir = "8nodes"
structures_path = joinpath(WRKDIR,"structures",structure_dir)
instances_path = joinpath(structures_path,"instances/.big_files")
instance_path_clust = joinpath(instances_path,instance_clust)
instance_ftr = "08n8760h_ftr"
instance_path_ftr = joinpath(instances_path,instance_ftr)
params_path_ftr = joinpath(instances_path,instance_ftr,".big_files","parameters")
results_dir = "results"

@info "Creating specifications"
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

@info "Create output directories path and directories"
if nosun && clust_method == "day" && specs.transmission
# Each case is specified below
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    results_clust_path = joinpath(instance_path_clust,string("output_tri_day_nosun"))
elseif nosun && clust_method == "day"
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    results_clust_path = joinpath(instance_path_clust,string("output_tri_nodal_day_nosun"))
elseif nosun && specs.transmission
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    results_clust_path = joinpath(instance_path_clust,string("output_tri_nosun"))
elseif nosun
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    results_clust_path = joinpath(instance_path_clust,string("output_tri_nodal_nosun"))
elseif clust_method == "day" && specs.transmission
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    results_clust_path = joinpath(instance_path_clust,string("output_tri_day"))
elseif clust_method == "day"
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    results_clust_path = joinpath(instance_path_clust,string("output_tri_nodal_day"))
elseif specs.transmission
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    results_clust_path = joinpath(instance_path_clust,string("output_tri"))
else
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    results_clust_path = joinpath(instance_path_clust,string("output_tri_nodal"))
end

output_path = results_clust_path * "_fix"
mkpath(output_path)
mkpath(joinpath(output_path,results_dir));

@info "Reading clust variables directories"
results_path = joinpath(results_clust_path,results_dir)
variables_clust = load(joinpath(results_path, "variables.jld2"));

@info "Reading parameters FTR"
# Argument nosun takes solar energy out of the optimisation"
# Updating the parameters to the clustering instance

## TODO: revise parameters (check if they are FTR)
parameters = change_time_parameters(params_path_ftr, instance_path_ftr)

@info "Creating the energy system model"
(model, VariablesDict, ObjectivesDict) = EnergySystemModel(parameters, specs)
# To write the model: write_to_file(model,joinpath(output_path, model_dir,"model.mps"))

@info "Adding capacity cuts"
# Fixing capacities
fix.(model[:p̄_gn], variables_clust["p̄_gn"]; force=true)
if specs.transmission
    fix.(model[:f̄_l], variables_clust["f̄_l"]; force=true)
end
if specs.storage
    fix.(model[:b̄_sn], variables_clust["b̄_sn"]; force=true)
end
if specs.hydro
    fix.(model[:h̄_hn], variables_clust["h̄_hn"]; force=true)
end

@info "Optimizing FTR model with capacity cuts"

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

@info "Saving results (JLD2)"
JLD2.save(joinpath(output_path, results_dir, "variables.jld2"), variables; compress = true)
JLD2.save(joinpath(output_path, results_dir, "objectives.jld2"), objectives; compress = true)
JLD2.save(joinpath(output_path, results_dir, "expressions.jld2"), expressions; compress = true)

@info string("Ended instance ",ARGS[1]," @ ",now())