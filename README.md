# EnergySystemModeling.jl
[![Docs Image](https://img.shields.io/badge/docs-latest-blue.svg)](https://gamma-opt.github.io/EnergySystemModeling.jl/dev/)
![Runtests](https://github.com/gamma-opt/EnergySystemModeling.jl/workflows/Runtests/badge.svg)

Julia library for solving the *generation and transmission capacity expansion problem*, implemented as a *linear program* using JuMP. On top of the optimisation model, the library provides:

- **Time series aggregation** (`src/aggreg.jl`) — chronological clustering of demand and renewable availability series, used to reduce the number of time periods in an instance while keeping the temporal order.
- **Plotting** (`src/plotting.jl`) — a set of standard result plots (dispatch, capacities, storage levels, transmission flows, loss of load).
- **Data generation** (`src/Data_generation/`) — building instance datasets for European regions from [GlobalEnergyGIS](https://github.com/niclasmattsson/GlobalEnergyGIS) output.

The documentation contains the full mathematical formulation of the model.

The library is authored by *Lucas Condeixa*, *Fabricio Oliveira*, and *Jaan Tollander de Balsch* in the Systems Analysis Laboratory at Aalto University.


## Model features
The constraints included in the model are selected through the `Specs` struct. All fields default to `false` except `transmission`.

| Field | Description |
|---|---|
| `transmission` | Transmission lines between nodes, with investment and flow limits |
| `renewable_target` | Minimum share of generation coming from renewable technologies |
| `carbon_cap` | Emission reduction target relative to a reference year |
| `nuclear_limit` | Maximum allowed share of nuclear generation |
| `storage` | Storage technologies (charge/discharge, capacity investment) |
| `ramping` | Ramp-up and ramp-down limits on generation |
| `voltage_angles` | DC power flow voltage angle constraints |
| `hydro` | Hydro reservoirs and run-of-river, including environmental flow |
| `hydro_simple` | Simplified hydro representation (use instead of `hydro`) |

Time periods may represent more than one hour (`τ_t`), which is what makes the aggregated instances produced by `src/aggreg.jl` usable directly in the model.


## Repository structure
```
src/                     Library source
  model.jl               Specs, Params, EnergySystemModel, Expressions
  io.jl                  Instance loading, JSON/JLD2 serialisation, clustered instances
  aggreg.jl              Chronological time series aggregation
  plotting.jl            Result plots
  Data_generation/       Instance generation from GlobalEnergyGIS data
examples/
  constants/             constants.json shared by all instances
  structures/            Node structures (8nodes, 11nodes, ...) and their instances
  run.jl, plot.jl        Example driver scripts
docs/                    Documenter.jl documentation sources
.triton/exe/             Slurm driver scripts used on the Triton cluster
```


## Installation
This library can be installed directly from GitHub
```
pkg> add https://github.com/gamma-opt/EnergySystemModeling.jl
```

Julia 1.10 or newer is required (see `[compat]` in `Project.toml`).


## Instance data
An instance is defined by two paths: the shared **constants** directory and the **instance** directory.

`constants/constants.json` holds the general parameters `kappa`, `mu`, `C`, `C_bar`, `C_E`, `R_E`, `r` and `Fmin`.

The instance directory contains:

- `indices.json` — fields `G`, `G_r`, `N`, `L`, `L_ind`, `T`, `S`, `H`
- `nodes/` — per-node time series (`1.csv`, `2.csv`, ...) with demand and availability columns
- `gen_technology.csv`, `gen_capacity.csv` — generation technology costs, lifetimes, efficiencies, emissions, ramp limits and initial capacities
- `transmission.csv` — line cost, distance, lifetime, capacity limits and susceptance
- `storage.csv`, `sto_capacity.csv` — storage technology parameters and capacity limits
- `hydro.csv`, `hydro_technology.csv`, `hydro_capacity.csv` — hydro reservoir and run-of-river parameters
- `nodes_specs.csv` — node-level specifications

Each instance directory carries a `README.md` listing the node, technology and line indexing for that instance (indices start at 1, following the Julia convention).


## Usage
`examples/run.jl` demonstrates the workflow on the instances shipped in `examples/structures`.

### Loading parameters and building the model
```julia
using EnergySystemModeling

constants_path = joinpath("examples", "constants")
instance_path = joinpath("examples", "structures", "8nodes", "instances", "ftr", "08n8760h_ftr")

parameters = Params(constants_path, instance_path)

specs = Specs(
    transmission=true,
    renewable_target=true,
    carbon_cap=true,
    nuclear_limit=false,
    storage=true,
    ramping=true,
    voltage_angles=false,
    hydro=true,
    hydro_simple=false
)

# EnergySystemModel returns the JuMP model together with dictionaries
# describing the variables and objective terms that were created.
(model, VariablesDict, ObjectivesDict) = EnergySystemModel(parameters, specs)
```

### Solving
```julia
using Gurobi, JuMP

optimizer = optimizer_with_attributes(Gurobi.Optimizer,
                                      "TimeLimit" => 60*60*2,
                                      "LogFile" => joinpath(output_path, "gurobi.log"))
set_optimizer(model, optimizer)
set_optimizer_attributes(model, "Method" => 2)
set_optimizer_attributes(model, "Crossover" => 0)
set_optimizer_attributes(model, "NumericFocus" => 1)

optimize!(model)
```

Barrier without crossover (`Method => 2`, `Crossover => 0`) is a good default for the large instances; `NumericFocus` and `ScaleFlag` help with the numerically harder ones.

### Extracting results
`JuMPVar` and `JuMPObj` turn the solved model into plain dictionaries of arrays, and `Expressions` computes derived quantities (renewable share `κ′`, hydro share `μ′`, emission reduction `C′_E`).

```julia
variables = JuMPVar(model, VariablesDict)
objectives = JuMPObj(model, ObjectivesDict)
expressions = Expressions(parameters, specs, variables)
```

### Saving and loading results
JSON, for small instances and for inspecting results by hand:
```julia
save_json(specs, joinpath(output_path, "specs.json"))
save_json(parameters, joinpath(output_path, "parameters.json"))
save_json(variables, joinpath(output_path, "variables.json"))
save_json(objectives, joinpath(output_path, "objectives.json"))
save_json(expressions, joinpath(output_path, "expressions.json"))

specs = load_json(Specs, joinpath(output_path, "specs.json"))
parameters = load_json(Params, joinpath(output_path, "parameters.json"))
```

JLD2, which is what the cluster runs use — considerably faster and smaller for full-year instances:
```julia
using JLD2
JLD2.save(joinpath(output_path, "variables.jld2"), variables; compress = true)
JLD2.save(joinpath(output_path, "objectives.jld2"), objectives; compress = true)
JLD2.save(joinpath(output_path, "expressions.jld2"), expressions; compress = true)
```

Parsing the CSV instance data is expensive, so it is worth doing once with `perform_Params`, which builds the `Params` object and stores it as `parameters.jld2`:
```julia
perform_Params(constants_path, instance_path, params_path)
```

### Plotting
Individual plotting functions (`plot_objective_values`, `plot_generation_dispatch`, `plot_generation_capacities`, `plot_generation_capacities_stacked`, `plot_transmission_flow`, `plot_transmission_capacities`, `plot_transmission_bars`, `plot_storage_level`, `plot_storage_capacities`, `plot_loss_of_load`, `plot_box`, `plot_box_all`, `plot_dispatch_bars`) each return a `Plots` object.

`perform_plotting` runs a selected subset of them and writes both PDF and PNG into `<plots_output_path>/pdf` and `<plots_output_path>/png`:
```julia
Plots_specs = Dict{String,Bool}(
    "p1" => true,    # objective function values
    "p2" => true,    # dispatch and storage levels, per node
    "p3" => true,    # storage capacities
    "p4" => true,    # generation dispatch levels (box plots)
    "p5" => true,    # generation capacities, stacked
    "p6" => false,   # consolidated dispatch vs demand
    "p7" => false,   # transmission flow, per line
    "p8" => false,   # transmission capacities
    "p9" => false,   # consolidated transmission flow
    "p10" => false,  # loss of load
    "p11" => false   # in development
)

perform_plotting(Plots_specs, parameters, variables, objectives, expressions, plots_output_path)
```
The `pdf` and `png` subdirectories must exist beforehand.


## Time series aggregation
Solving a full 8760-hour instance is expensive, so the library provides chronological aggregation that merges neighbouring time steps into clusters while preserving their order. Each resulting period keeps a weight, which the model reads as `τ_t`.

The aggregation runs on a matrix `series` whose columns are the individual time series (demand per node, followed by availability per technology and node) and whose rows are time steps:

```julia
# Instance carrying the series and the aggregation settings
_SeriesInstance = load_series_instance(
    series,
    block_size,        # clusters merged per iteration
    current_k,         # starting number of clusters
    stopping_k,        # target number of clusters
    dm,                # discrepancy metric, e.g. :ed (euclidean) or :wd (ward)
    rep_value          # cluster representative value, e.g. :mean
)

# Initial clustering state (one cluster per time step)
_ClustInstance = load_clust_instance(k_cent, series_clust)

# Aggregate down to stopping_k, recording every intermediate clustering
find_clusters!(_SeriesInstance, _ClustInstance, _DistUpdate)
```

`write_clust_instance!` then materialises the clustered instances as instance directories that `Params`/`change_time_parameters` can read back, and `read_clust_instance` / `read_clusters` load them again.

To solve an aggregated instance, load the pre-computed full-resolution parameters and swap in the clustered time dimension:
```julia
parameters = change_time_parameters(params_path_ftr, instance_path_clust)
# optionally: nosun = true to zero out solar availability
```
`break_FTR_parameters_in_periods` does the same for a sub-range of periods, which is used for the representative-day runs.

`.triton/exe/aggregation/` contains the driver scripts for this workflow: `form_aggreg_instances_new.jl` runs the aggregation itself, and `generate_clust_instances.jl` writes out the resulting instance folders.


## Data generation
`src/Data_generation/` builds instance datasets for European countries and regions from [GlobalEnergyGIS](https://github.com/niclasmattsson/GlobalEnergyGIS) output. `run_data_generation.jl` is the entry point — set the SSP scenario, target year, ERA year, region list, time horizon and number of technologies at the top of the file, then call `create_data_sets`, which writes a new instance directory under `examples/structures`. See `docs/src/datageneration.md` for details.


## Running on a cluster
`.triton/exe/` holds the Slurm scripts used to run the model on Aalto's Triton cluster, split into `aggregation/` (building clustered instances), `opt/` (solving them) and `plot/` (post-processing). The `opt` directory contains variants of the driver script for the different experiment setups: `run.jl`, `run_FTR.jl` (full time resolution), `run_clust.jl` (aggregated instances), `run_FTR_days.jl` (representative days), `run_fix.jl` and `run_min_cap.jl`. These scripts contain absolute paths and are meant to be adapted rather than run as-is.


## Development
Install the [Julia](https://julialang.org/) programming language.

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

Gurobi is a powerful commercial optimizer that provides a free academic license. We can interface with Gurobi in Julia using [`Gurobi.jl`](https://github.com/jump-dev/Gurobi.jl). Here are the steps to install Julia and Gurobi to run the program:

1) Obtain a license of *Gurobi* and install Gurobi solver by following the instructions on [Gurobi's website](http://www.gurobi.com/).

2) Make sure the `GUROBI_HOME` environmental variable is set to the path of the Gurobi directory. This is part of standard installation. The Gurobi library will be searched for in `GUROBI_HOME/lib` on Unix platforms and `GUROBI_HOME\bin` on Windows. If the library is not found, check that your version is listed in `deps/build.jl`. The environmental variable can be set by appending `export GUROBI_HOME="<path>/gurobi811/linux64"` to `.bashrc` file. Replace the `<path>`, platform `linux64` and version number `811` with the values of your Gurobi installation.

3) Install `Gurobi.jl` in Julia's package manager by running commands
   ```
   pkg> add Gurobi
   pkg> build Gurobi
   ```


## Documentation
The project documentation is created using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/). It covers the mathematical formulation (`index.md`), the aggregation method (`aggreg.md`), plotting (`plotting.md`), data generation (`datageneration.md`) and the full API reference (`api.md`).

To build the documentation, navigate inside the `docs` directory and run the command
```bash
julia make.jl
```
