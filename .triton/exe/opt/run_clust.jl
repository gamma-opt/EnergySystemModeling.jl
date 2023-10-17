## Code for running optimisation linked it to the clustering technique
using Logging
using EnergySystemModeling
using JLD2
using JuMP
using Gurobi
using Dates
using DelimitedFiles

# Define optimisation time limit
opt_tlim = 60*60*72

@info string("Starting instance ",ARGS[1]," @ ",now())

# Instance including clustering info
instance_clust = ARGS[1]                # Clustering instance with info
clust_method = ARGS[2]                  # Clustering method (e.g., "","dc"/"dc_new"/"dc_nosun"/"day_RMSE")
nosun = parse(Bool,ARGS[3])             # Boolean regarding solar availability

# NOTE: Change the output folder's name accordingly to the instance
instance_dir = string("output_tri_",clust_method)
results_dir = "results"

WRKDIR = "/scratch/work/condeil1/EnergySystemModeling.jl/examples"
cd(WRKDIR)

@info "Defining constants and instance files path"
constants_dir = "constants"
constants_path = joinpath(WRKDIR,constants_dir)
structure_dir = "8nodes"
structures_path = joinpath(WRKDIR,"structures",structure_dir)
instances_path = joinpath(structures_path,"instances/.big_files")
instance_path_clust = joinpath(instances_path,instance_clust)
instance_ftr = "08n8760h_ftr"
params_path_ftr = joinpath(instances_path,instance_ftr,".big_files","parameters");

@info "Clustering parameters"
aggregation_path = "/scratch/work/condeil1/EnergySystemModeling.jl/.triton/exe/aggregation/"

# Defining the clustering results path
if string(instance_clust[end-1]) == "m"
        if occursin("day",clust_method)
                clust_method_path = joinpath(aggregation_path,string("aggreg_output_ed_",clust_method),string("clust_out_ed_",clust_method,".jld2"));
        else
                clust_method_path = joinpath(aggregation_path,string("aggreg_output_ed",clust_method),string("clust_out_ed",clust_method,".jld2"));
        end
elseif string(instance_clust[end-1]) == "w"
        if occursin("day",clust_method)
                clust_method_path = joinpath(aggregation_path,string("aggreg_output_wd_",clust_method),string("clust_out_wd_",clust_method,".jld2"));
        else
                clust_method_path = joinpath(aggregation_path,string("aggreg_output_wd",clust_method),string("clust_out_wd",clust_method,".jld2"));
        end
end

@info "Creating specifications"
# Model specifications
# NOTE: to run nodal, disable transmission
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

@info "Loading parameters to instance"
# Updating the parameters to the clustering instance
ParamsDict = load(joinpath(params_path_ftr,"parameters.jld2"));
parameters = read_clust_instance(clust_method, clust_method_path, instance_clust, ParamsDict ; nosun=nosun);

if occursin("day",clust_method)
        # Number of representative hours from the instance_clust
        rep_hours = parse(Int,instance_clust[9:12])

        # Number of days from the instance_clust
        hours_day = 24
        days_year = 365

        @assert rep_hours%hours_day == 0 "Non-integer number of days"
        num_days = rep_hours/hours_day |> Int
        @assert num_days < days_year "More days than possible in an year"
end

@info "Creating output directories"
if nosun && occursin("day",clust_method) && specs.transmission
# Each case is specified below
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    output_path = joinpath(instance_path_clust,string("output_tri_day_nosun"))
elseif nosun && occursin("day",clust_method)
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    output_path = joinpath(instance_path_clust,string("output_tri_nodal_day_nosun"))
elseif nosun && specs.transmission
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    output_path = joinpath(instance_path_clust,string("output_tri_nosun"))
elseif nosun
    # No solar series, rep. days method, connected network, standard clustering method (Ward's)
    output_path = joinpath(instance_path_clust,string("output_tri_nodal_nosun"))
elseif occursin("day",clust_method) && occursin("RMSE",clust_method) && specs.transmission
    # Solar series included, rep. days method, connected network, standard clustering method (Ward's)
    output_path = joinpath(instance_path_clust,string("output_tri_day_RMSE"))
elseif occursin("day",clust_method) && specs.transmission
    # Solar series included, rep. days method, connected network, standard clustering method (Ward's)
    output_path = joinpath(instance_path_clust,string("output_tri_day"))
elseif occursin("day",clust_method)
    # Solar series included, rep. days method, connected network, standard clustering method (Ward's)
    output_path = joinpath(instance_path_clust,string("output_tri_nodal_day"))
elseif specs.transmission
    # Solar series included, rep. days method, connected network, standard clustering method (Ward's)
    output_path = joinpath(instance_path_clust,string("output_tri"))
else
    # Solar series included, rep. days method, connected network, standard clustering method (Ward's)
    output_path = joinpath(instance_path_clust,string("output_tri_nodal"))
end

@info string("Output directory path created/updated: ", output_path)
mkpath(output_path)
results_dir = "results"
mkpath(joinpath(output_path,results_dir));

@info "Creating energy system model"
(model, VariablesDict, ObjectivesDict) = EnergySystemModel(parameters, specs)
# write_to_file(model,joinpath(output_path, model_dir,"model.mps"))

if occursin("day",clust_method)
        # Adjust storage levels for rep. days method
        @info "Adjusting storage levels and hydro reservoirs levels to rep. days"
        T_sto = sort(vcat(((collect(2:num_days+1).-1)*hours_day),((collect(1:num_days-1)*hours_day).-1)))

        @constraint(model,
            s0_rd[s in parameters.S, n in parameters.N, t in T_sto],
            model[:b_snt][s,n,t]/1000 == 0.5*model[:b̄_sn][s,n]/1000)

        @constraint(model, 
            h0_rd[h in parameters.H, n in parameters.N, t in T_sto],
            model[:w_hnt][h,n,t]/1000 == parameters.Wmax_hn[h,n]/1000)
end;

@info "Optimizing the model"
optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => opt_tlim,
                                      "LogFile" => joinpath(output_path, "gurobi.log"))
set_optimizer(model, optimizer)
set_optimizer_attributes(model, "NumericFocus" => 1)
set_optimizer_attributes(model, "ScaleFlag" => 2)

# Scaling: ("ScaleFlag" => [-1,3])
# Numeric focus: ("NumericFocus" => [0,3])
# set_optimizer_attributes(model, "Method" => 0)
# set_optimizer_attributes(model, "Crossover" => 1)

optimize!(model);

@info "Extracting results"
variables = JuMPVar(model, VariablesDict)
objectives = JuMPObj(model, ObjectivesDict)
expressions = Expressions(parameters, specs, variables);

@info "Saving results (JLD2)"
JLD2.save(joinpath(output_path, results_dir, "variables.jld2"), variables; compress = true)
JLD2.save(joinpath(output_path, results_dir, "objectives.jld2"), objectives; compress = true)
JLD2.save(joinpath(output_path, results_dir, "expressions.jld2"), expressions; compress = true);

@info string("Ended instance ",ARGS[1]," @ ",now());