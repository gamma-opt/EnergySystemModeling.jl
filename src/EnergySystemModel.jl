module EnergySystemModel

using JuMP, JSON, CSV, DataFrames

export Specs, energy_system_model, load_parameters

"""Load parameters"""
function load_parameters(instance_path)
    # Load indexes and constant parameters
    indices = JSON.parsefile(joinpath(instance_path, "indices.json"))

    # Indices
    G = indices["G"]
    G_r = indices["G_r"]
    N = indices["N"]
    L = indices["L"]
    # FIXME: clustering, time steps, τ parameter
    T = 1:indices["T"]
    S = indices["S"]

    # Constant parameters
    constants = JSON.parsefile(joinpath(instance_path, "constants.json"))
    κ = constants["kappa"]
    C = constants["C"]

    # Load time clustered parameters
    D = zeros(length(T), length(N))
    A_wind = zeros(length(T), length(N))
    A_solar = zeros(length(T), length(N))
    for n in N
        # Load node values from CSV files.
        df = CSV.read(joinpath(instance_path, "nodes", "$n.csv")) |> DataFrame
        # TODO: clustering
        D[:, n] = df.Dem_Inc[1] .* df.Load_mod .* df.Max_Load[1]
        A_wind[:, n] = df.Avail_Win
        A_solar[:, n] = df.Avail_Sol
    end

    # Load technology parameters
    technology = joinpath(instance_path, "technology.csv") |>
        CSV.read |> DataFrame
    I_g = technology.I
    M_g = technology.M
    C_g = technology.C
    r⁻ = technology.r_minus
    r⁺ = technology.r_plus

    # Load transmission parameters
    transmission = joinpath(instance_path, "transmission.csv") |>
        CSV.read |> DataFrame
    I_l = transmission.I
    M_l = transmission.M
    C_l = transmission.C
    B_l = transmission.B

    # Load storage parameters
    storage = joinpath(instance_path, "storage.csv") |>
        CSV.read |> DataFrame
    ξ_s = storage.xi
    I_s = storage.I
    C_s = storage.C
    b⁰_sn = storage[:, [Symbol("b0_$n") for n in N]]

    # Return tuple (maybe namedtuple?)
    return (G, G_r, N, L, T, S, κ, C, D, A_wind, A_solar, I_g, M_g, C_g, r⁻,
            r⁺, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b⁰_sn)
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
