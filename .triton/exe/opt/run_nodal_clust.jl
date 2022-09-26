## Code for running the nodal optimisation linking it to the clustering technique
using Logging
using EnergySystemModeling
using JLD2
using JuMP
using Gurobi
using Dates
using DelimitedFiles

@info string("Starting instance ",ARGS[1]," @ ",now())

# Instance including clustering info
instance_clust = ARGS[1]                # Clustering instance with info
clust_method = ARGS[2]                  # Clustering method ("","dc"/"dc_new"/"dc_nosun")
nosun = parse(Bool,ARGS[3])             # Boolean regarding solar availability

WRKDIR = "/scratch/work/condeil1/EnergySystemModeling.jl/examples"
cd(WRKDIR)

# Instance as a string(number)
instance = parse(Int,instance_clust[9:12]) |> string

@info "Loading parameters"
constants_dir = "constants"
constants_path = joinpath(WRKDIR,constants_dir)
structure_dir = "8nodes"
structures_path = joinpath(WRKDIR,"structures",structure_dir)
instances_path = joinpath(structures_path,"instances/.big_files")
instance_path_clust = joinpath(instances_path,instance_clust)
instance_ftr = "08n8760h_ftr"
params_path_ftr = joinpath(instances_path,instance_ftr,".big_files","parameters")

@info "Clustering parameters"
aggregation_path = "/scratch/work/condeil1/EnergySystemModeling.jl/.triton/exe/aggregation/"
# Defining the clustering results path
if string(instance_clust[end-1]) == "m"
        clust_method_path = joinpath(aggregation_path,string("aggreg_output_ed",clust_method),string("clust_out_ed",clust_method,".jld2"));
elseif string(instance_clust[end-1]) == "w"
        clust_method_path = joinpath(aggregation_path,string("aggreg_output_wd",clust_method),string("clust_out_wd",clust_method,".jld2"));
end

@info "Creating specifications"
# Model specifications
specs = Specs(
        transmission=false,
        renewable_target=true,
        carbon_cap=true,
        nuclear_limit=false,
        storage=true,
        ramping=true,
        voltage_angles=false,
        hydro=true,
        hydro_simple=false
)

@info "Adjusting parameters to instance"
# Updating the parameters to the clustering instance
ParamsDict = load(joinpath(params_path_ftr,"parameters.jld2"));
parameters = read_clust_instance(clust_method_path, instance, ParamsDict ; nosun=nosun);

@info "Creating output directories"
if nosun
        if clust_method == ""
                output_path = joinpath(instance_path_clust,string("output_tri_nodal_nosun"))
        else
                output_path = joinpath(instance_path_clust,string("output_tri_nodal_",clust_method))
        end
else
        if clust_method == ""
                output_path = joinpath(instance_path_clust,string("output_tri_nodal"))
        else
                output_path = joinpath(instance_path_clust,string("output_tri_nodal_",clust_method))
        end
end
mkpath(output_path)
results_dir = "results"
mkpath(joinpath(output_path,results_dir))

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

@info "Saving results (JLD2)"
JLD2.save(joinpath(output_path, results_dir, "variables.jld2"), variables; compress = true)
JLD2.save(joinpath(output_path, results_dir, "objectives.jld2"), objectives; compress = true)
JLD2.save(joinpath(output_path, results_dir, "expressions.jld2"), expressions; compress = true)

@info string("Ended instance ",ARGS[1]," @ ",now())