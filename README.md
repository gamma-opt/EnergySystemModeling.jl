# EnergySystemModeling.jl
[![Docs Image](https://img.shields.io/badge/docs-latest-blue.svg)](https://gamma-opt.github.io/EnergySystemModeling.jl/dev/)
![Runtests](https://github.com/gamma-opt/EnergySystemModeling.jl/workflows/Runtests/badge.svg)

Julia library for solving the *transmission capacity expansion problem*, implemented as *linear program* using JuMP. The documentation contains more details about the model.

The library is authored by *Lucas Condeixa*, *Fabricio Oliveira*, and *Jaan Tollander de Balsch* in Systems Analysis Laboratory in Aalto university.


## Usage
Inside the `examples` directory, we have [`run.jl`](./examples/run.jl) file, which demonstrates the usage of this library by running the example [instance](./examples/instance).

```julia
using Gurobi, JuMP
using EnergySystemModeling

# Load parameters.
parameters = Params(joinpath("examples", "instance"))

# Define specs.
specs = Specs(
    renewable_target=true,
    storage=true,
    ramping=false,
    voltage_angles=false
)

# Create the model.
model = EnergySystemModel(parameters, specs)

# Optimizer the model using Gurobi optimizer.
optimizer = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "TimeLimit" => 5*60
)
optimize!(model, optimizer)

# Extract values from the model.
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

We recommend to check out the documentation for plotting.


## Installation
This library can be installed directly from GitHub
```
pkg> add https://github.com/gamma-opt/EnergySystemModeling.jl
```


## Development
Install [Julia](https://julialang.org/) programming language.

Clone the repository
```bash
git clone https://github.com/gamma-opt/EnergySystemModeling.jl.git
```

In the project root directory, install packages locally using Julia's package manager.
```
pkg> dev .
```

Install a solver such as Gurobi.


## Installing Solver
It's up to the user to choose a suitable solver for solving the JuMP model. For small instances, GLPK is sufficient, but for large instances, we recommend commercial solvers such as Gurobi or CPLEX.

Gurobi is a powerful commercial optimizer that provides a free academic license. We can interface with Gurobi in Julia using [`Gurobi.jl`](https://github.com/JuliaOpt/Gurobi.jl). Here are the steps to install Julia and Gurobi to run the program:

1) Obtain a license of *Gurobi* and install Gurobi solver by following the instructions on [Gurobi's website](http://www.gurobi.com/).

2) Make sure the `GUROBI_HOME` environmental variable is set to the path of the Gurobi directory. This is part of standard installation. The Gurobi library will be searched for in `GUROBI_HOME/lib` on Unix platforms and `GUROBI_HOME\bin` on Windows. If the library is not found, check that your version is listed in `deps/build.jl`. The environmental variable can be set by appending `export GUROBI_HOME="<path>/gurobi811/linux64"` to `.bashrc` file. Replace the `<path>`, platform `linux64` and version number `811` with the values of your Gurobi installation.

3) Install `Gurobi.jl` in Julia's package manager by running commands
   ```
   pkg> add Gurobi
   pkg> build Gurobi
   ```


## Documentation
The project documentation is created using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/). To build the documentation, navigate inside the `docs` directory and run the command
```bash
julia make.jl
```
