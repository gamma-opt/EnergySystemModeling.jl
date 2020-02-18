using Logging, JuMP

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModel

@info "Loading parameters"
parameters = load_parameters("instance")
specs = Specs(true, true, false, false)

@info "Creating the energy system model"
model = energy_system_model(parameters, specs)

@info "Optimizing the model"
# using GLPK
# optimizer = with_optimizer(GLPK.Optimizer, tm_lim = 60000)
using Gurobi
optimizer = with_optimizer(Gurobi.Optimizer, TimeLimit=5*60)
optimize!(model, optimizer)
