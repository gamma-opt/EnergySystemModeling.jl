using Logging, JuMP

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModeling

@info "Create output directory"
output = "output"
if !ispath(output)
    mkdir(output)
end

@info "Loading parameters"
parameters = Parameters("instance")
specs = Specs(
    renewable_target=true,
    storage=true,
    ramping=false,
    voltage_angles=false
)

@info "Creating the energy system model"
model = energy_system_model(parameters, specs)

@info "Optimizing the model"
using Gurobi
optimizer = with_optimizer(Gurobi.Optimizer, TimeLimit=5*60,
                           LogFile=joinpath(output, "gurobi.log"))
optimize!(model, optimizer)

@info "Extract results"
variables = Variables(model)
objectives = Objectives(model)

@info "Save results"
save_results(specs, parameters, variables, objectives, output)

@info "Plotting"
using Plots

savefig(plot_objective_values(objectives),
        joinpath(output, "objectives.svg"))

for n in parameters.N
    savefig(plot_generation_dispatch(parameters, variables, n),
            joinpath(output, "generation_dispatch_n$n.svg"))
    savefig(plot_generation_capacities(parameters, variables, n),
            joinpath(output, "generation_capacities_n$n.svg"))
    savefig(plot_storage(parameters, variables, n),
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
