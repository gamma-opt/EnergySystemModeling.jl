using Logging

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModeling

@info "Create output directory"
output = "output"
mkpath(output)

@info "Loading parameters"
parameters = Params(joinpath("examples","Smallinstance"))
specs = Specs(
    renewable_target=true,
    carbon_cap=true,
    nuclear_limit=false,
    storage=true,
    ramping=true,
    voltage_angles=false
)

@info "Creating the energy system model"
model = EnergySystemModel(parameters, specs)

@info "Optimizing the model"
using Gurobi, JuMP
optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60*60*2,
                                      "LogFile" => joinpath(output, "gurobi.log"))
set_optimizer(model, optimizer)
set_optimizer_attributes(model, "Method" => 2)
set_optimizer_attributes(model, "Crossover" => 0)
set_optimizer_attributes(model, "NumericFocus" => 1)
optimize!(model)

@info "Extract results"
variables = Variables(model)
objectives = Objectives(model)
expressions = Expressions(parameters, variables)

@info "Save results"
save_json(specs, joinpath(output, "specs.json"))
save_json(parameters, joinpath(output, "parameters.json"))
save_json(variables, joinpath(output, "variables.json"))
save_json(objectives, joinpath(output, "objectives.json"))

@info "Plotting"
using Plots
using StatsPlots
pyplot()

savefig(plot_objective_values(objectives),
        joinpath(output, "objectives.pdf"))

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
savefig(plot_storage_capacities(parameters, variables, expressions),
        joinpath(output, "storage_capacities.pdf"))

savefig(plot_box_all(parameters, variables, expressions),
        joinpath(output, "boxplotall.pdf"))

savefig(plot_dispatch_bars(parameters, variables, expressions),
        joinpath(output, "dispatchbars.pdf"))

for l in 1:length(parameters.L)
    savefig(plot_transmission_flow(parameters, variables, expressions, l),
            joinpath(output, "transmission_flow_l$l.pdf"))
end

savefig(plot_transmission_capacities(parameters, variables, expressions),
        joinpath(output, "transmission_capacities.pdf"))

savefig(plot_transmission_bars(parameters, variables, expressions),
        joinpath(output, "transmission_bars.pdf"))

savefig(plot_loss_of_load(parameters, variables, expressions),
        joinpath(output, "loss_of_load.pdf"))

getdispatch(joinpath(output))