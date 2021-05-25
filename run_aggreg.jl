using Logging

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModeling

@info "Creating aggreg_TS directory"
output_dir = "examples//aggreg_out"

mkpath(output_dir)

@info "Loading parameters"
constants_path = "examples//constants"
structure = "8nodes"
structures_path = joinpath("examples//structures",structure)
instance = "small"
instances_path = joinpath(structures_path,"instances",instance)

parameters = Params(constants_path, instances_path)

parameters.D_nt
