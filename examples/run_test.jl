## pwd()
## cd("C://...//EnergySystemModel.jl//examples")
##
using Logging

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModeling

## Creating directory
@info "Create output directory"
output = "output_test"
mkpath(output)

## Loading parameters
@info "Loading parameters"
parameters = Params("smallinstance")
specs = Specs(
    renewable_target=true,
    carbon_cap=true,
    nuclear_limit=false,
    storage=true,
    ramping=true,
    voltage_angles=false
)

## Creating model object
@info "Creating the energy system model"
model = EnergySystemModel(parameters, specs)

## Optimizing
@info "Optimizing the model"
using Gurobi, JuMP
optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60*60*2,
                                      "LogFile" => joinpath(output, "gurobi.log"))
set_optimizer(model, optimizer)
set_optimizer_attributes(model, "Method" => 2)
set_optimizer_attributes(model, "Crossover" => 0)
set_optimizer_attributes(model, "NumericFocus" => 1)
optimize!(model)

## Storing results
@info "Extract results"
variables = Variables(model)
objectives = Objectives(model)
expressions = Expressions(parameters, variables)

## Saving results
@info "Save results"
save_json(specs, joinpath(output, "specs.json"))
save_json(parameters, joinpath(output, "parameters.json"))
save_json(variables, joinpath(output, "variables.json"))
save_json(objectives, joinpath(output, "objectives.json"))

## Plotting
@info "Plotting"
using Plots
using StatsPlots
pyplot()

## Plotting part 1
@info "Plotting OF"
savefig(plot_objective_values(objectives),
        joinpath(output, "objectives.pdf")
)

## Plotting part 2
@info "Plotting dispatch and storage levels"

for n in parameters.N
    savefig(plot_generation_dispatch(parameters, variables, expressions, n),
            joinpath(output, "generation_dispatch_n$n.pdf"))
    savefig(plot_generation_capacities(parameters, variables, expressions, n),
            joinpath(output, "generation_capacities_n$n.pdf"))
    savefig(plot_storage_level(parameters, variables, expressions, n),
            joinpath(output, "storage_n$n.pdf"))
    savefig(plot_box(parameters, variables, expressions, n),
           joinpath(output, "boxplot$n.pdf"))
    
end

## Plotting part 3
@info "Plotting storage capacities"
savefig(plot_storage_capacities(parameters, variables, expressions),
        joinpath(output, "storage_capacities.pdf")
)

## Plotting part 4
@info "Plotting dispatch all"
savefig(plot_box_all(parameters, variables, expressions),
        joinpath(output, "boxplotall.pdf")
)

## Plotting part 5
@info "Plotting dispatch all (boxplot)"
savefig(plot_dispatch_bars(parameters, variables, expressions),
        joinpath(output, "dispatchbars.pdf")
)

## Plotting part 6
@info "Plotting transmission flow"
for l in 1:length(parameters.L)
    savefig(plot_transmission_flow(parameters, variables, expressions, l),
            joinpath(output, "transmission_flow_l$l.pdf"))
end

## Plotting part 7
@info "Plotting transmission capacities"
savefig(plot_transmission_capacities(parameters, variables, expressions),
        joinpath(output, "transmission_capacities.pdf")
)

## Plotting part 8
@info "Plotting transmission flow"
savefig(plot_transmission_bars(parameters, variables, expressions),
        joinpath(output, "transmission_bars.pdf")
)

## Plotting part 9
@info "Plotting LoL"
savefig(plot_loss_of_load(parameters, variables, expressions),
        joinpath(output, "loss_of_load.pdf")
)

## Plotting in development
@info "Plotting stacked dispatch"