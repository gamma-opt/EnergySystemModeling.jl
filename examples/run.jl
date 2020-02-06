using Logging, JuMP

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModel

@info "Loading parameters"
parameters = load_parameters("instance")
specs = Specs(true, true, true, true)

@info "Creating the energy system model"
model = energy_system_model(parameters, specs)

# @info "Optimizing the model"
# using GLPK
# optimizer = with_optimizer(GLPK.Optimizer)
# optimize!(model, optimizer)
