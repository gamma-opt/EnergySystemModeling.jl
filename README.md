# EnergySystemModeling.jl
Julia library for solving the *transmission capacity expansion problem*, implemented as *linear program* using JuMP.

The library is authored by *Lucas Condeixa*, *Fabricio Oliveira*, and *Jaan Tollander de Balsch* in Systems Analysis Laboratory in Aalto university.


## Usage
Inside the `examples` directory, we have [`run.jl`](./examples/run.jl) file, which demonstrates the usage of this library by running the example [instance](./examples/instance).

```julia
using EnergySystemModeling

parameters = Params(joinpath("examples", "instance"))
specs = Specs(
    renewable_target=true,
    storage=true,
    ramping=false,
    voltage_angles=false
)
model = EnergySystemModel(parameters, specs)

using Gurobi, JuMP
optimizer = with_optimizer(Gurobi.Optimizer, TimeLimit=5*60)
optimize!(model, optimizer)

variables = Variables(model)
objectives = Objectives(model)
```

Saving values to JSON.
```julia
save_json(specs, joinpath("output", "specs.json"))
save_json(parameters, joinpath("output", "parameters.json"))
save_json(variables, joinpath("output", "variables.json"))
save_json(objectives, joinpath("output", "objectives.json"))
```

Loading values from JSON.
```julia
specs = load_json(Specs, joinpath("output", "specs.json"))
parameters = load_json(Params, joinpath("output", "parameters.json"))
variables = load_json(Variables, joinpath("output", "variables.json"))
objectives = load_json(Objectives, joinpath("output", "objectives.json"))
```

We can generate plots using the following functions.
```julia
plot_objective_values(objectives)
for n in parameters.N
    plot_generation_dispatch(parameters, variables, n)
    plot_generation_capacities(parameters, variables, n)
    plot_storage_level(parameters, variables, n)
    plot_storage_capacities(parameters, variables, n)
end
for l in 1:length(parameters.L)
    plot_transmission_flow(parameters, variables, l)
end
plot_transmission_capacities(parameters, variables)
plot_loss_of_load(parameters, variables)
```


## Installation
This library can be installed directly from GitHub
```
pkg> add https://github.com/jaantollander/EnergySystemModeling.jl
```


## Development
Install [Julia](https://julialang.org/) programming language.

In the project root directory, install packages locally using Julia's package manager.
```
pkg> activate .
pkg> instantiate
```

Before we can run the `examples/run.jl` script, we require an optimizer for solving the linear program. One such solver is [GLPK](https://github.com/JuliaOpt/GLPK.jl), an open-source solver for linear programming. We can install it using Julia's package manager.
```
pkg> add GLPK
```

Alternatively, you could choose to install a different solver such as Gurobi, a powerful commercial solver. In this case, you need to modify the code inside `run.jl` file. Now, we can run the script.
```bash
julia run.jl
```
