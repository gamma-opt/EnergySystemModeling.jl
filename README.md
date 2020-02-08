# EnergySystemModel.jl
Julia library for solving the *transmission capacity expansion problem*, implemented as *linear program* using JuMP.

The library is authored by *Lucas Condeixa*, *Fabricio Oliveira*, and *Jaan Tollander de Balsch* in Systems Analysis Laboratory in Aalto university.


## Usage
Inside the `examples` directory, we have [`run.jl`](./examples/run.jl) file, which demonstrates the usage of this library by running the example [instance](./examples/instance).


## Installation
This library can be installed directly from GitHub
```
pkg> add https://github.com/jaantollander/EnergySystemModel.jl
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

Alternatively, you could choose to install a different solver. In this case, you need to modify the code inside `run.jl` file. Now, we can run the script.
```bash
julia run.jl
```
