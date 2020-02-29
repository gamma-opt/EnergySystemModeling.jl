# API
## Model
```@docs
EnergySystemModel
Specs
Params
Variables
Objectives
EnergySystemModel(::Params, ::Specs)
```

## Input
```@docs
equivalent_annual_cost
Params(::AbstractString)
```

## Output
```@docs
Variables(::EnergySystemModel)
Objectives(::EnergySystemModel)
save_results
```

## Plotting
```@docs
plot_objective_values
plot_generation_dispatch
plot_generation_capacities
plot_transmission_flow
plot_transmission_capacities
plot_storage
plot_storage_capacities
plot_loss_of_load
```
