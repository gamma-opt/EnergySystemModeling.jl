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

## IO
```@docs
equivalent_annual_cost
Params(::AbstractString)
Variables(::EnergySystemModel)
Objectives(::EnergySystemModel)
save_json
load_json
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
