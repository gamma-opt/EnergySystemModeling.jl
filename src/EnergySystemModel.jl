module EnergySystemModel

using JuMP, JSON, CSV, DataFrames

export Specs, energy_system_model, load_parameters

"""Load parameters"""
function load_parameters(instance_path)
    # Load indexes and constant parameters
    indices = JSON.parsefile(joinpath(instance_path, "indices.json"))

    # Indices. Convert JSON values to right types.
    G = indices["G"] |> Array{Int}
    G_r = indices["G_r"] |> Array{Int}
    N = indices["N"] |> Array{Int}
    # TODO: return L_pairs
    L_pairs = indices["L"] |> Array{Array{Int}}
    L = 1:length(L_pairs)
    # TODO: clustering, time steps, τ parameter
    τ = 1
    T = 1:indices["T"]
    S = indices["S"] |> Array{Int}

    # Constant parameters
    constants = JSON.parsefile(joinpath(instance_path, "constants.json"))
    κ = constants["kappa"]
    C = constants["C"]

    # Load time clustered parameters
    D_nt = zeros(length(N), length(T))
    A1_nt = zeros(length(N), length(T)) # wind
    A2_nt = zeros(length(N), length(T)) # solar
    for n in N
        # Load node values from CSV files.
        df = CSV.read(joinpath(instance_path, "nodes", "$n.csv")) |> DataFrame
        # TODO: clustering
        D_nt[n, :] = df.Dem_Inc[1] .* df.Load_mod .* df.Max_Load[1]
        A1_nt[n, :] = df.Avail_Win
        A2_nt[n, :] = df.Avail_Sol
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
    # TODO: convert to array
    b⁰_sn = storage[:, [Symbol("b0_$n") for n in N]]

    # Return tuple (maybe namedtuple?)
    return (G, G_r, N, L, T, S, κ, C, τ, A1_nt, A2_nt, D_nt, I_g, M_g,
            C_g, r⁻, r⁺, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b⁰_sn)
end

"""Specs"""
struct Specs
    # TODO: booleans
end

"""Create energy system model."""
function energy_system_model(
            G, G_r, N, L, T, S, κ, C, τ, A1_nt, A2_nt, D_nt, I_g, M_g,
            C_g, r⁻, r⁺, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b⁰_sn)::Model
    # Create an instance of JuMP model.
    model = Model()

    # -- Variables --
    @variable(model, p_gnt[g in G, n in N, t in T]≥0)
    @variable(model, p̄_gn[g in G, n in N]≥0)
    @variable(model, σ_nt[n in N, t in T]≥0)
    @variable(model, f_lt[l in L, t in T])
    @variable(model, f̄_l[l in L])
    @variable(model, b_snt[s in S, n in N, t in T]≥0)
    @variable(model, b̄_sn[s in S, n in N]≥0)
    @variable(model, b⁺_snt[s in S, n in N, t in T]≥0)
    @variable(model, b⁻_snt[s in S, n in N, t in T]≥0)
    @variable(model, θ_nt[n in N, t in T]≥0)
    @variable(model, θ′_nt[n in N, t in T]≥0)

    # -- Objective --
    @expression(model, f1,
        sum((I_g[g]+M_g[g])*p̄_gn[g,n] for g in G, n in N))
    @expression(model, f2,
        sum(C_g[g]*p_gnt[g,n,t]*τ for g in G, n in N, t in T))
    @expression(model, f3,
        sum(C*σ_nt[n,t]*τ for n in N, t in T))
    @expression(model, f4,
        sum((I_l[l]+M_l[l])*f̄_l[l] for l in L))
    # TODO: linerize abs
    @expression(model, f5,
        sum(C_l[l]*f_lt[l,t]*τ for l in L, t in T))
    @expression(model, f6,
        sum(I_s[s]*b̄_sn[s,n] for s in S, n in N))
    @expression(model, f7,
        sum(C_s[s]*(b⁺_snt[s,n,t] + b⁻_snt[s,n,t])*τ for s in S, n in N, t in T))
    @objective(model, Min, f1 + f2 + f3 + f4 + f5 + f6 + f7)

    # -- Constraints --


    return model
end

end # module
