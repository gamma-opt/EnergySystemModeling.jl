using JuMP
using Gurobi

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModel

# TODO: parameter values

model = energy_system_model()

optimize!(model, with_optimizer(Gurobi.Optimizer))
