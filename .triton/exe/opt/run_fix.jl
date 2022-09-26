using Logging
using EnergySystemModeling
using JLD2
using JuMP
using Gurobi
using Dates
using DelimitedFiles

@info string("Starting instance (fixed capacity) ",ARGS[1]," @ ",now()) 

WRKDIR = "/scratch/work/condeil1/EnergySystemModeling.jl/examples"
cd(WRKDIR)

@info "Loading instance specs"
constants_dir = "constants"
constants_path = joinpath(WRKDIR,constants_dir)
structure_dir = "8nodes"
structures_path = joinpath(WRKDIR,"structures",structure_dir)
instances_path = joinpath(structures_path,"instances/.big_files")
instance_clust = ARGS[1]
instance_path_clust = joinpath(instances_path,instance_clust)
instance_ftr = "08n8760h_ftr"
instance_path_ftr = joinpath(instances_path,instance_ftr)

output_ftr_path = joinpath(instance_path_ftr,"output_tri")
output_clust_path = joinpath(instance_path_clust,"output_tri")
output_fix_path = joinpath(instance_path_clust,"output_fix")

results_dir = "results"

mkpath(output_fix_path)
mkpath(joinpath(output_fix_path,results_dir))

params_path_ftr = joinpath(instances_path,instance_ftr,".big_files","parameters")

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

@info "Reading parameters (FTR)"
ParamsDict = load(joinpath(params_path_ftr, "parameters.jld2"));
parameters_ftr = Params(
    ParamsDict["region_n"], ParamsDict["max_dem_n"], ParamsDict["technology_g"], ParamsDict["G"], ParamsDict["G_r"], ParamsDict["N"], ParamsDict["L"], ParamsDict["L_ind"], ParamsDict["T"], 
    ParamsDict["S"], ParamsDict["H"], ParamsDict["κ"], ParamsDict["μ"], ParamsDict["C"], ParamsDict["C̄"], ParamsDict["C_E"], ParamsDict["R_E"], ParamsDict["τ_t"], ParamsDict["Gmin_gn"], ParamsDict["Gmax_gn"],
    ParamsDict["A_gnt"], ParamsDict["D_nt"], ParamsDict["I_g"], ParamsDict["M_g"], ParamsDict["C_g"], ParamsDict["e_g"], ParamsDict["E_g"], ParamsDict["r⁻_g"],
    ParamsDict["r⁺_g"], ParamsDict["I_l"], ParamsDict["M_l"], ParamsDict["C_l"], ParamsDict["B_l"], ParamsDict["e_l"], ParamsDict["Tmin_l"], ParamsDict["Tmax_l"],
    ParamsDict["ξ_s"], ParamsDict["I_s"], ParamsDict["C_s"], ParamsDict["Smin_sn"], ParamsDict["Smax_sn"], ParamsDict["Wmax_hn"], ParamsDict["Wmin_hn"], ParamsDict["Hmax_hn"],
    ParamsDict["Hmin_hn"], ParamsDict["HRmax_n"], ParamsDict["Fmin_n"], ParamsDict["AH_nt"], ParamsDict["AR_nt"], ParamsDict["I_h"], ParamsDict["M_h"], ParamsDict["C_h"],
    ParamsDict["e_h"], ParamsDict["E_h"], ParamsDict["r⁻_h"], ParamsDict["r⁺_h"]
)

@info "Extracting clust. variables (clust)"
variables = load(joinpath(output_clust_path, results_dir, "variables.jld2"));
# variables = load(joinpath(output_ftr_path, results_dir, "variables.jld2"));

@info "Creating the FTR model (fixed capacity)"
(model, VariablesDict, ObjectivesDict) = EnergySystemModel(parameters_ftr, specs)

optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60*60*48,
                                    "LogFile" => joinpath(output_fix_path, "gurobi.log"))
set_optimizer(model, optimizer)

# Aggressive scaling: ("ScaleFlag" => 2)
set_optimizer_attributes(model, "NumericFocus" => 1)
set_optimizer_attributes(model, "ScaleFlag" => 2)

## Other possible attributes:
# set_optimizer_attributes(model, "Crossover" => 0)
# set_optimizer_attributes(model, "Method" => 2)
# set_optimizer_attributes(model, "NumericFocus" => 1)

@info "Adding capacity cuts"
# Fixing capacities
fix.(model[:p̄_gn], variables["p̄_gn"]; force=true)
if specs.transmission
    fix.(model[:f̄_l], variables["f̄_l"]; force=true)
end
if specs.storage
    fix.(model[:b̄_sn], variables["b̄_sn"]; force=true)
end
if specs.hydro
    fix.(model[:h̄_hn], variables["h̄_hn"]; force=true)
end

@info "Optimizing FTR model with capacity cuts"
optimize!(model)

@info "Extracting FTR results"
variables_fix = JuMPVar(model, VariablesDict)
objectives_fix = JuMPObj(model, ObjectivesDict)
expressions_fix = Expressions(parameters_ftr, specs, variables_fix)

# @info "Saving results (JSON)"
# save_json(specs, joinpath(output_fix_path, results_dir, "specs.json"))
# save_json(variables_fix, joinpath(output_fix_path, results_dir, "variables.json"))
# save_json(objectives_fix, joinpath(output_fix_path, results_dir, "objectives.json"))
# save_json(expressions_fix, joinpath(output_fix_path, results_dir, "expressions.json"))

@info "Saving results (JLD2)"
JLD2.save(joinpath(output_fix_path, results_dir, "variables.jld2"), variables_fix; compress = true)
JLD2.save(joinpath(output_fix_path, results_dir, "objectives.jld2"), objectives_fix; compress = true)
JLD2.save(joinpath(output_fix_path, results_dir, "expressions.jld2"), expressions_fix; compress = true)

## Write a out file
# open(joinpath(instance_path_clust,"out2.txt"), "w") do io
#     writedlm(io, "")
# end

@info string("Ended instance (fixed capacity) ",ARGS[1]," @ ",now()) 