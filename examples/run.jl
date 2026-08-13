using Logging
push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModeling

ENV["GKSwstype"]="nul"                      # Prevent opening plots windows (must be set before Plots loads)
cd(@__DIR__)

@info "Loading parameters"
constants_path = "constants"
structure = "8nodes"
structures_path = joinpath("structures",structure)
# Instance directory, relative to structures/<structure>/instances.
# Defaults to the 168-hour instance, which solves in a couple of minutes; pass
# another as the first argument, e.g. `julia --project=. examples/run.jl ftr/08n8760h_ftr`.
instance = isempty(ARGS) ? joinpath(".big_files","08n0168h_ftr") : ARGS[1]
instances_path = joinpath(structures_path,"instances",instance)
isdir(instances_path) || error("Instance directory not found: $(abspath(instances_path))")

@info "Creating output directory"
output_dir = joinpath(instances_path,"output_local")
results_dir = "results"
plots_dir = "plots"
csv_dir = "csv"

mkpath(output_dir)
mkpath(joinpath(output_dir,results_dir))
mkpath(joinpath(output_dir,plots_dir))
mkpath(joinpath(output_dir,plots_dir,"pdf"))
mkpath(joinpath(output_dir,plots_dir,"png"))
mkpath(joinpath(output_dir,csv_dir))

parameters = Params(constants_path, instances_path)
specs = Specs(
        transmission=true,
        renewable_target=true,
        carbon_cap=true,
        nuclear_limit=false,
        storage=true,
        ramping=true,
        voltage_angles=false,
        hydro=true
)

@info "Creating the energy system model"
(model, VariablesDict, ObjectivesDict) = EnergySystemModel(parameters, specs)

@info "Optimizing the model"
using Gurobi, JuMP
optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60*60*2,
                                      "LogFile" => joinpath(output_dir, "gurobi.log"))
set_optimizer(model, optimizer)
set_optimizer_attributes(model, "Method" => 2)
set_optimizer_attributes(model, "NumericFocus" => 1)
# Crossover is disabled in the cluster runs (see .triton/exe/opt) to save time on the
# full-year instances. Locally it is worth keeping: without it the barrier hits
# "Numerical trouble encountered" on these instances and returns no solution.
# set_optimizer_attributes(model, "Crossover" => 0)

optimize!(model)

if !has_values(model)
    error("No solution available: termination status $(termination_status(model)). " *
          "See $(joinpath(output_dir, "gurobi.log")).")
end

@info "Extracting results"
variables = JuMPVar(model, VariablesDict)
objectives = JuMPObj(model, ObjectivesDict)
expressions = Expressions(parameters, specs, variables)

@info "Saving results (JSON)"
save_json(specs, joinpath(output_dir, results_dir, "specs.json"))
save_json(parameters, joinpath(output_dir, results_dir, "parameters.json"))
save_json(variables, joinpath(output_dir, results_dir, "variables.json"))
save_json(objectives, joinpath(output_dir, results_dir, "objectives.json"))

@info "Plotting"
using Plots
using StatsPlots
gr()

## Plotting part 1: Objective function values
@info "Plotting OF"
savefig(plot_objective_values(objectives, parameters.τ_t, parameters.T),
        joinpath(output_dir, plots_dir, "pdf", "objectives.pdf"))
savefig(plot_objective_values(objectives, parameters.τ_t, parameters.T),
        joinpath(output_dir, plots_dir, "png", "objectives.png"))

## Plotting part 2: Dispatch and storage
@info "Plotting dispatch and storage levels"

for n in parameters.N
    savefig(plot_generation_dispatch(parameters, variables, expressions, n),
            joinpath(output_dir, plots_dir,"pdf","generation_dispatch_n$n.pdf"))
    savefig(plot_generation_capacities(parameters, variables, expressions, n),
            joinpath(output_dir, plots_dir,"pdf","generation_capacities_n$n.pdf"))
    savefig(plot_storage_level(parameters, variables, expressions, n),
            joinpath(output_dir, plots_dir,"pdf","storage_n$n.pdf"))
    savefig(plot_box(parameters, variables, expressions, n),
           joinpath(output_dir, plots_dir,"pdf","boxplot$n.pdf"))
    savefig(plot_generation_dispatch(parameters, variables, expressions, n),
           joinpath(output_dir, plots_dir,"png","generation_dispatch_n$n.png"))
    savefig(plot_generation_capacities(parameters, variables, expressions, n),
           joinpath(output_dir, plots_dir,"png","generation_capacities_n$n.png"))
    savefig(plot_storage_level(parameters, variables, expressions, n),
           joinpath(output_dir, plots_dir,"png","storage_n$n.png"))
    savefig(plot_box(parameters, variables, expressions, n),
          joinpath(output_dir, plots_dir,"png","boxplot$n.png"))

end

## Plotting part 3: Storage capacity
@info "Plotting storage capacities"

savefig(plot_storage_capacities(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"pdf","storage_capacities.pdf")
        )
savefig(plot_storage_capacities(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"png","storage_capacities.png")
        )

## Plotting part 4: Generation technologies dispatch levels
@info "Plotting dispatch all"

savefig(plot_box_all(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"pdf","boxplotall.pdf"))
savefig(plot_box_all(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"png","boxplotall.png"))

## Plotting part 5: Generation capacities (including hydro)
@info "Plotting capacities stacked"

savefig(plot_generation_capacities_stacked(parameters, variables, expressions),
            joinpath(output_dir, plots_dir,"pdf","generation_capacities_stacked.pdf"))
        savefig(plot_generation_capacities_stacked(parameters, variables, expressions),
            joinpath(output_dir, plots_dir,"png","generation_capacities_stacked.png"))

## Plotting part 6: Consolidated dispatch vs demand
@info "Plotting dispatch all (boxplot)"

savefig(plot_dispatch_bars(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"pdf","dispatchbars.pdf"))
savefig(plot_dispatch_bars(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"png","dispatchbars.png"))

## Plotting part 7: Transmission flow (per line)
@info "Plotting transmission flow"

for l in 1:length(parameters.L_ind)
    savefig(plot_transmission_flow(parameters, variables, expressions, l),
            joinpath(output_dir, plots_dir,"pdf","transmission_flow_L$l.pdf"))
    savefig(plot_transmission_flow(parameters, variables, expressions, l),
            joinpath(output_dir, plots_dir,"png","transmission_flow_L$l.png"))

end

## Plotting part 8: Transmission capacities
@info "Plotting transmission capacities"

savefig(plot_transmission_capacities(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"pdf","transmission_capacities.pdf"))
savefig(plot_transmission_capacities(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"png","transmission_capacities.png"))

## Plotting part 9: Consolidated transmission flow
@info "Plotting transmission flow"

savefig(plot_transmission_bars(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"pdf","transmission_bars.pdf"))
savefig(plot_transmission_bars(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"png","transmission_bars.png"))

## Plotting part 10: Lost of load
@info "Plotting LoL"

savefig(plot_loss_of_load(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"pdf","loss_of_load.pdf"))
savefig(plot_loss_of_load(parameters, variables, expressions),
        joinpath(output_dir, plots_dir,"png","loss_of_load.png"))

## Plotting in development
@info "Plotting stacked dispatch (not ready)"