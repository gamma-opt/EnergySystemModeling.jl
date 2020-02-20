using Logging, JuMP

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModel

@info "Create output directory"
output = "output"
if !ispath(output)
    mkdir(output)
end

@info "Loading parameters"
parameters = load_parameters("instance")
specs = Specs(true, true, false, false)

@info "Creating the energy system model"
model = energy_system_model(parameters, specs)

@info "Optimizing the model"
# using GLPK
# optimizer = with_optimizer(GLPK.Optimizer, tm_lim = 60000)
using Gurobi
optimizer = with_optimizer(Gurobi.Optimizer, TimeLimit=5*60,
                           LogFile=joinpath(output, "gurobi.log"))
optimize!(model, optimizer)

@info "Plotting"
using Plots

savefig(plot_objective_values(model::Model),
        joinpath(output, "objectives.svg"))

for n in parameters.N
    savefig(plot_generation_dispatch(parameters, model, n),
            joinpath(output, "generation_dispatch_n$n.svg"))
    savefig(plot_generation_capacities(parameters, model, n),
            joinpath(output, "generation_capacities_n$n.svg"))
    savefig(plot_storage(parameters, model, n),
            joinpath(output, "storage_n$n.svg"))
    savefig(plot_storage_capacities(parameters, model, n),
           joinpath(output, "storage_capacities_n$n.svg"))
end

for l in 1:length(parameters.L)
    savefig(plot_transmission_flow(parameters, model, l),
            joinpath(output, "transmission_flow_l$l.svg"))
end

savefig(plot_transmission_capacities(parameters, model),
        joinpath(output, "transmission_capacities.svg"))

savefig(plot_loss_of_load(parameters, model),
        joinpath(output, "loss_of_load.svg"))
