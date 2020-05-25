using Logging

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModeling

@info "Create output directory"
output = "output"
mkpath(output)

@info "Loading parameters"
parameters = Params("instance")
specs = Specs(
    renewable_target=true,
    storage=true,
    ramping=false,
    voltage_angles=false
)

@info "Creating the energy system model"
model = EnergySystemModel(parameters, specs)

@info "Optimizing the model"
using Gurobi, JuMP
optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 10*60,
                                      "LogFile" => joinpath(output, "gurobi.log"))
set_optimizer(model, optimizer)
optimize!(model)

@info "Extract results"
variables = Variables(model)
objectives = Objectives(model)

@info "Save results"
save_json(specs, joinpath(output, "specs.json"))
save_json(parameters, joinpath(output, "parameters.json"))
save_json(variables, joinpath(output, "variables.json"))
save_json(objectives, joinpath(output, "objectives.json"))

@info "Plotting"
using Plots

savefig(plot_objective_values(objectives),
        joinpath(output, "objectives.svg"))

for n in parameters.N
    savefig(plot_generation_dispatch(parameters, variables, n),
            joinpath(output, "generation_dispatch_n$n.svg"))
    savefig(plot_generation_capacities(parameters, variables, n),
            joinpath(output, "generation_capacities_n$n.svg"))
    savefig(plot_storage_level(parameters, variables, n),
            joinpath(output, "storage_n$n.svg"))
    savefig(plot_storage_capacities(parameters, variables, n),
           joinpath(output, "storage_capacities_n$n.svg"))
end

for l in 1:length(parameters.L)
    savefig(plot_transmission_flow(parameters, variables, l),
            joinpath(output, "transmission_flow_l$l.svg"))
end

savefig(plot_transmission_capacities(parameters, variables),
        joinpath(output, "transmission_capacities.svg"))

savefig(plot_loss_of_load(parameters, variables),
        joinpath(output, "loss_of_load.svg"))
