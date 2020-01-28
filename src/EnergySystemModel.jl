module EnergySystemModel

using JuMP, JSON, CSV
using JuMP.Containers

export Specs, energy_system_model, load_parameters

"""Load parameters"""
function load_parameters(instance_path)
    # Load indexes and constant parameters
    p = JSON.parsefile(joinpath(instance_path, "parameters.json"))

    constants = p["constants"]
    indices = p["indices"]
    G = indices["G"]
    G_r = indices["G_r"]
    N = indices["N"]
    L = indices["L"]
    S = indices["S"]

    # TODO: time steps
    T = []

    # Time clustered parameters
    # We use DenseAxisArray because indices might be non-integers.
    demand = DenseAxisArray(zeros(length(T), length(N)), T, N)
    wind = DenseAxisArray(zeros(length(T), length(N)), T, N)
    solar = DenseAxisArray(zeros(length(T), length(N)), T, N)
    for n in N
        df = CSV.read(joinpath(instance_path, "$n.csv"))
        D = Float64.(df.Load_mod)
        W = Float64.(df.Avail_Win)
        S = Float64.(df.Avail_Sol)
        max_load = Float64.(df.Max_Load[1])
        dem_inc = Float64.(df.Dem_Inc[1])

        # TODO: clustering?

        # Compute demands
        demand[:, n] = dem_inc .* D[T] .* max_load
        wind[:, n] = W[T]
        solar[:, n] = S[T]
    end

    # Load technology parameters
    technology = CSV.read(joinpath(instance_path, "technology.csv"))

    # Load transmission parameters
    transmission = CSV.read(joinpath(instance_path, "transmission.csv"))

    # Load storage parameters
    storage = CSV.read(joinpath(instance_path, "storage.csv"))

    # Return namedtuple
    return ()
end

"""Specs"""
struct Specs
    # TODO: booleans
end

"""Create energy system model."""
function energy_system_model()::Model
    # Create an instance of JuMP model.
    model = Model()

    # TODO: variables
    # TODO: objective
    # TODO: constraints

    return model
end

end # module
