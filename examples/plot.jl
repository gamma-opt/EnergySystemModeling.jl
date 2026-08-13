using Logging
push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModeling, JSON

ENV["GKSwstype"]="nul"                      # Prevent opening plots windows (must be set before Plots loads)
cd(@__DIR__)

# Re-creates the plots from the results saved by run.jl, without re-solving the model.
# Takes the same instance argument as run.jl, e.g.
# `julia --project=. examples/plot.jl ftr/08n8760h_ftr`.
@info "Loading results"
structure = "8nodes"
structures_path = joinpath("structures",structure)
instance = isempty(ARGS) ? joinpath(".big_files","08n0168h_ftr") : ARGS[1]
instances_path = joinpath(structures_path,"instances",instance)

output_dir = joinpath(instances_path,"output_local")
results_path = joinpath(output_dir,"results")
plots_path = joinpath(output_dir,"plots")
isdir(results_path) || error("No results found in $(abspath(results_path)). Run examples/run.jl first.")

mkpath(joinpath(plots_path,"pdf"))
mkpath(joinpath(plots_path,"png"))

"""Rebuild an array from the nested arrays written by `save_json`. JSON nests the
outermost dimension last, so e.g. `p_gnt` is stored as `[t][n][g]` and has to be
stacked back into `(g,n,t)`."""
function nested_array(x::AbstractVector)
    isempty(x) && return Float64[]
    first(x) isa AbstractVector || return Float64.(x)
    stack(nested_array(xi) for xi in x)
end

specs = load_json(Specs, joinpath(results_path, "specs.json"))
parameters = load_json(Params, joinpath(results_path, "parameters.json"))
variables = Dict{String, Array{Float64}}(
    k => nested_array(v) for (k, v) in JSON.parsefile(joinpath(results_path, "variables.json")))
objectives = Dict{String, Float64}(JSON.parsefile(joinpath(results_path, "objectives.json")))

# Expressions are derived from the results rather than stored, so they are recomputed here.
expressions = Expressions(parameters, specs, variables)

## Plotting specifications
Plots_specs = Dict{String,Bool}(
    "p1" => true,     # objective function values
    "p2" => true,     # dispatch and storage levels (per node)
    "p3" => true,     # storage capacities
    "p4" => true,     # generation dispatch levels (box plots)
    "p5" => true,     # generation capacities, stacked
    "p6" => true,     # consolidated dispatch vs demand
    "p7" => true,     # transmission flow (per line)
    "p8" => true,     # transmission capacities
    "p9" => true,     # consolidated transmission flow
    "p10" => true,    # loss of load
    "p11" => false    # in development
)

perform_plotting(Plots_specs, parameters, variables, objectives, expressions, plots_path)

@info "Plots written to $(abspath(plots_path))"
