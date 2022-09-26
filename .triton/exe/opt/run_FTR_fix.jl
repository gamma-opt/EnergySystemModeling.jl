using Logging
using EnergySystemModeling
using JLD2
using JuMP
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
instances_path = joinpath(structures_path,"instances")
instance_clust = ARGS[1]
instance_path_clust = joinpath(instances_path,instance_clust)
instance_ftr = "08n8760h_ftr"

@info "Creating output directories"
output_fix_path = joinpath(instance_path_clust,"output_fix_FTR")
mkpath(output_fix_path)
results_dir = "results"

params_path_ftr = joinpath(instances_path,instance_ftr,".big_files","parameters")
output_FTR_path = joinpath(instances_path,instance_ftr,"output_tri")


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

@info "Extracting FTR variables"
variables_ftr = load(joinpath(output_FTR_path, results_dir, "variables.jld2"));

@info "Creating the energy system model"
(model, VariablesDict, ObjectivesDict) = EnergySystemModel(parameters, specs)

@info "Optimizing the model"
using Gurobi
optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60*60*48,
                                      "LogFile" => joinpath(output_fix_path, "gurobi.log"))
set_optimizer(model, optimizer)
set_optimizer_attributes(model, "NumericFocus" => 1)
set_optimizer_attributes(model, "ScaleFlag" => 2)

# Scaling: ("ScaleFlag" => [-1,3])
# Numeric focus: ("NumericFocus" => [0,3])
# set_optimizer_attributes(model, "Method" => 0)
# set_optimizer_attributes(model, "Crossover" => 0)

@info "Adding capacity cuts"
fix.(model[:p̄_gn], variables_ftr["p̄_gn"]; force=true)
if specs.transmission
    fix.(model[:f̄_l], variables_ftr["f̄_l"]; force=true)
end
if specs.storage
    fix.(model[:b̄_sn], variables_ftr["b̄_sn"]; force=true)
end
if specs.hydro
    fix.(model[:h̄_hn], variables_ftr["h̄_hn"]; force=true)
end

optimize!(model)

@info "Extracting results"
variables = JuMPVar(model, VariablesDict)
objectives = JuMPObj(model, ObjectivesDict)
expressions = Expressions(parameters, specs, variables)

@info "Saving results (JLD2)"
JLD2.save(joinpath(output_fix_path, results_dir, "variables.jld2"), variables; compress = true)
JLD2.save(joinpath(output_fix_path, results_dir, "objectives.jld2"), objectives; compress = true)
JLD2.save(joinpath(output_fix_path, results_dir, "expressions.jld2"), expressions; compress = true)

open(joinpath(instance_path_clust,"out_fix.txt"), "w") do io
        writedlm(io, "")
end

@info string("Ended instance ",ARGS[1]," @ ",now())