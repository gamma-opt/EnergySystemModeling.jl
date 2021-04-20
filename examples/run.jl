using Logging

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModeling

@info "Creating output directory"
output_dir = "output_small"
results_dir = "results"
plots_dir = "plots"
csv_dir = "csv"

mkpath(output_dir)
mkpath(joinpath(output_dir,results_dir))
mkpath(joinpath(output_dir,plots_dir))
mkpath(joinpath(output_dir,csv_dir))

@info "Loading parameters"
constants_path = "constants"
structure = "8nodes"
structures_path = joinpath("structures",structure)
instance = "small"
instances_path = joinpath(structures_path,"instances",instance)

parameters = Params(constants_path, instances_path)
specs = Specs(
    renewable_target=true,
    carbon_cap=true,
    nuclear_limit=false,
    storage=true,
    ramping=true,
    voltage_angles=false,
    hydro = true
)

@info "Creating the energy system model"
model = EnergySystemModel(parameters, specs)

@info "Optimizing the model"
using Gurobi, JuMP
optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60*60*2,
                                      "LogFile" => joinpath(output_dir, "gurobi.log"))
set_optimizer(model, optimizer)
set_optimizer_attributes(model, "Method" => 2)
set_optimizer_attributes(model, "Crossover" => 0)
set_optimizer_attributes(model, "NumericFocus" => 1)
optimize!(model)

@info "Extracting results"
variables = Variables(model)
objectives = Objectives(model)
expressions = Expressions(parameters, variables)

@info "Saving results (JSON)"
save_json(specs, joinpath(output_dir, results_dir, "specs.json"))
save_json(parameters, joinpath(output_dir, results_dir, "parameters.json"))
save_json(variables, joinpath(output_dir, results_dir, "variables.json"))
save_json(objectives, joinpath(output_dir, results_dir, "objectives.json"))

@info "Plotting"
using Plots
ENV["GKSwstype"]="nul"
using StatsPlots
gr()

## Plotting part 1
@info "Plotting OF"
savefig(plot_objective_values(objectives),
        joinpath(output_dir, plots_dir, "objectives.pdf")
        )

## Plotting part 2
@info "Plotting dispatch and storage levels"

for n in parameters.N
    savefig(plot_generation_dispatch(parameters, variables, expressions, n),
            joinpath(output_dir, plots_dir, "generation_dispatch_n$n.pdf"))
    savefig(plot_generation_capacities(parameters, variables, expressions, n),
            joinpath(output_dir, plots_dir, "generation_capacities_n$n.pdf"))
    savefig(plot_storage_level(parameters, variables, expressions, n),
            joinpath(output_dir, plots_dir, "storage_n$n.pdf"))
    savefig(plot_box(parameters, variables, expressions, n),
           joinpath(output_dir, plots_dir, "boxplot$n.pdf"))
end

## Plotting part 3
@info "Plotting storage capacities"

savefig(plot_storage_capacities(parameters, variables, expressions),
        joinpath(output_dir, plots_dir, "storage_capacities.pdf")
        )

## Plotting part 4
@info "Plotting dispatch all"

savefig(plot_box_all(parameters, variables, expressions),
        joinpath(output_dir, plots_dir, "boxplotall.pdf"))

## Plotting part 5
@info "Plotting capacities stacked"

savefig(plot_generation_capacities_stacked(parameters, variables, expressions),
            joinpath(output_dir, plots_dir, "generation_capacities_stacked.pdf"))

## Plotting part 6
@info "Plotting dispatch all (boxplot)"

savefig(plot_dispatch_bars(parameters, variables, expressions),
        joinpath(output_dir, plots_dir, "dispatchbars.pdf"))

## Plotting part 7
@info "Plotting transmission flow"

for l in 1:length(parameters.L_ind)
    savefig(plot_transmission_flow(parameters, variables, expressions, l),
            joinpath(output_dir, plots_dir, "transmission_flow_L$l.pdf"))
end

## Plotting part 8
@info "Plotting transmission capacities"

savefig(plot_transmission_capacities(parameters, variables, expressions),
        joinpath(output_dir, plots_dir, "transmission_capacities.pdf"))

## Plotting part 9
@info "Plotting transmission flow"

savefig(plot_transmission_bars(parameters, variables, expressions),
        joinpath(output_dir, plots_dir, "transmission_bars.pdf"))

## Plotting part 10
@info "Plotting LoL"

savefig(plot_loss_of_load(parameters, variables, expressions),
        joinpath(output_dir, plots_dir, "loss_of_load.pdf"))

## Plotting in development
@info "Plotting stacked dispatch (not ready)"