push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModeling, Parameters, JSON

output = "output"

parameters = load_json(EnergySystemModeling.Params, joinpath(output, "parameters.json"))
variables = load_json(EnergySystemModeling.Variables, joinpath(output, "variables.json"))
objectives = load_json(EnergySystemModeling.Objectives, joinpath(output, "objectives.json"))
expressions = Expressions(parameters, variables)

using Plots
using StatsPlots
pyplot()

savefig(plot_objective_values(objectives),
        joinpath(output, "objectives.svg"))

for n in parameters.N
    savefig(plot_generation_dispatch(parameters, variables, expressions, n),
            joinpath(output, "generation_dispatch_n$n.svg"))
    savefig(plot_generation_capacities(parameters, variables, expressions, n),
            joinpath(output, "generation_capacities_n$n.svg"))
    savefig(plot_storage_level(parameters, variables, expressions, n),
            joinpath(output, "storage_n$n.svg"))
    savefig(plot_box(parameters, variables, expressions, n),
           joinpath(output, "boxplot$n.svg"))
    
end
savefig(plot_storage_capacities(parameters, variables, expressions),
        joinpath(output, "storage_capacities.svg"))

savefig(plot_box_all(parameters, variables, expressions),
        joinpath(output, "boxplotall.svg"))

savefig(plot_dispatch_bars(parameters, variables, expressions),
        joinpath(output, "dispatchbars.svg"))

for l in 1:length(parameters.L)
    savefig(plot_transmission_flow(parameters, variables, expressions, l),
            joinpath(output, "transmission_flow_l$l.svg"))
end

savefig(plot_transmission_capacities(parameters, variables, expressions),
        joinpath(output, "transmission_capacities.svg"))

savefig(plot_transmission_bars(parameters, variables, expressions),
        joinpath(output, "transmission_bars.svg"))

savefig(plot_loss_of_load(parameters, variables, expressions),
        joinpath(output, "loss_of_load.svg"))