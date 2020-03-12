# API
## Model
```@docs
EnergySystemModel
Specs
Params
Variables
Objectives
Variables(::EnergySystemModel)
Objectives(::EnergySystemModel)
EnergySystemModel(::Params, ::Specs)
```

## IO
!!! note
    JSON stores multi-dimensional arrays as nested arrays. If we load an array `a` whose elements are accessed `a[i, j, k]`, the elements of the nested array are accessed in reverse order `a[k][j][i]`. However, the function [`load_json`](@ref) converts the nested arrays back to multi-dimensional arrays.

We use simple plain text formats CSV and JSON for input and output of numerical values.

```@docs
equivalent_annual_cost
Params(::AbstractString)
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
plot_storage_level
plot_storage_capacities
plot_loss_of_load
```
