# EnergySystemModel.jl
Julia library for solving the *transmission capacity expansion problem*, implemented as optimization model using JuMP.

The library is authored by *Lucas Condeixa*, *Fabricio Oliveira*, and *Jaan Tollander de Balsch* in Systems Analysis Laboratory in Aalto university.


## Installation
This library can be installed directly from GitHub
```
pkg> add https://github.com/jaantollander/EnergySystemModel.jl
```


## Development
Install [Julia](https://julialang.org/) programming language.

In the project root directory, install packages locally using
```
pkg> activate .
pkg> instantiate
```

In `examples` diretory there is `run.jl` script for running the [example instance](./examples/instance). However, to run the instance, we require an optimizer, such as [GLPK](https://github.com/JuliaOpt/GLPK.jl). Install your preferred optimizer, modify the script to use this optmizer and then run it using:
```bash
julia run.jl
```
