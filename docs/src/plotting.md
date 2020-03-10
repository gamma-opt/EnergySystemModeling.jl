# Plotting
```@example plots
using EnergySystemModeling
specs = load_json(Specs, joinpath("output", "specs.json"))
parameters = load_json(Params, joinpath("output", "parameters.json"))
variables = load_json(Variables, joinpath("output", "variables.json"))
objectives = load_json(Objectives, joinpath("output", "objectives.json"))
nothing; #hide
```

## Objectives
```@example plots
plot_objective_values(objectives)
```

## Generation Dispatch
```@example plots
n = 1
plot_generation_dispatch(parameters, variables, n)
```

## Generation Capacities
```@example plots
plot_generation_capacities(parameters, variables, n)
```

## Loss of Load
```@example plots
plot_loss_of_load(parameters, variables)
```

## Transmission Flow
```@example plots
l = 1
plot_transmission_flow(parameters, variables, l)
```

## Transmission Capacities
```@example plots
plot_transmission_capacities(parameters, variables)
```

## Storage Level
```@example plots
plot_storage_level(parameters, variables, n)
```

## Storage Capacities
```@example plots
plot_storage_capacities(parameters, variables, n)
```
