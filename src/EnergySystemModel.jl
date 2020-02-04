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
    L = indices["L"] |> Array{Array{Int}}
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
    r⁻_g = technology.r_minus
    r⁺_g = technology.r_plus

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
            C_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b⁰_sn)
end

"""Specs"""
struct Specs
    # TODO: booleans
end

"""Create energy system model."""
function energy_system_model(
            G, G_r, N, L, T, S, κ, C, τ, A1_nt, A2_nt, D_nt, I_g, M_g,
            C_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b⁰_sn)::Model
    # Create an instance of JuMP model.
    model = Model()

    # Indices of lines L
    L′ = 1:length(L)

    # TODO: CS_P / τ parameter

    ## -- Variables --
    @variable(model, p_gnt[g in G, n in N, t in T]≥0)
    @variable(model, p̄_gn[g in G, n in N]≥0)
    @variable(model, σ_nt[n in N, t in T]≥0)
    @variable(model, f_lt[l in L′, t in T])
    @variable(model, f_lt_abs[l in L′, t in T])
    @variable(model, f̄_l[l in L′])
    @variable(model, b_snt[s in S, n in N, t in T]≥0)
    @variable(model, b̄_sn[s in S, n in N]≥0)
    @variable(model, b⁺_snt[s in S, n in N, t in T]≥0)
    @variable(model, b⁻_snt[s in S, n in N, t in T]≥0)
    @variable(model, θ_nt[n in N, t in T]≥0)
    @variable(model, θ′_nt[n in N, t in T]≥0)

    ## -- Objective --
    @expression(model, f1,
        sum((I_g[g]+M_g[g])*p̄_gn[g,n] for g in G, n in N))
    @expression(model, f2,
        sum(C_g[g]*p_gnt[g,n,t]*τ for g in G, n in N, t in T))
    @expression(model, f3,
        sum(C*σ_nt[n,t]*τ for n in N, t in T))
    @expression(model, f4,
        sum((I_l[l]+M_l[l])*f̄_l[l] for l in L′))
    @expression(model, f5,
        sum(C_l[l]*f_lt_abs[l,t]*τ for l in L′, t in T))
    @expression(model, f6,
        sum(I_s[s]*b̄_sn[s,n] for s in S, n in N))
    @expression(model, f7,
        sum(C_s[s]*(b⁺_snt[s,n,t] + b⁻_snt[s,n,t])*τ for s in S, n in N, t in T))
    @objective(model, Min, f1 + f2 + f3 + f4 + f5 + f6 + f7)

    ## -- Constraints --
    # Energy balance (t=1)
    @constraint(model,
        [s in S, n in N, t in [1]],
        sum(p_gnt[g,n,t] for g in G) +
        σ_nt[n,t] +
        sum(f_lt[l,t] for (l,l′) in enumerate(L) if l′[1]==n) -
        sum(f_lt[l,t] for (l,l′) in enumerate(L) if l′[2]==n) +
        ξ_s[s] * b_snt[s,n,t] ==
        D_nt[n,t])

    # Energy balance (t>1)
    @constraint(model,
        [s in S, n in N, t in T[T.>1]],
        sum(p_gnt[g,n,t] for g in G) +
        σ_nt[n,t] +
        sum(f_lt[l,t] for (l,l′) in enumerate(L) if l′[1]==n) -
        sum(f_lt[l,t] for (l,l′) in enumerate(L) if l′[2]==n) +
        ξ_s[s] * (b_snt[s,n,t] - b_snt[s,n,t-1]) ==
        D_nt[n,t])

    # TODO: Generation capacity, GCN, parameter A

    # Minimum renewables share
    @constraint(model,
        sum(p_gnt[g,n,t] for g in G_r, n in N, t in T) ≥
        κ * sum(p_gnt[g,n,t] for g in G, n in N, t in T))

    # Shedding upper bound
    @constraint(model,
        [n in N, t in T],
        σ_nt[n,t] ≤ C * D_nt[n,t])

    # Transmission capacity
    @constraint(model,
        [l in L′, t in T],
        f_lt[l,t]≤f̄_l[l])
    @constraint(model,
        [l in L′, t in T],
        f_lt[l,t]≥-f̄_l[l])

    # Absolute value of transmission
    @constraint(model,
        [l in L′, t in T],
        f_lt_abs[l,t]≥f_lt[l,t])
    @constraint(model,
        [l in L′, t in T],
        f_lt_abs[l,t]≥-f_lt[l,t])

    # Charge and discharge (t=1)
    @constraint(model,
        [s in S, n in N, t in [1]],
        b⁺_snt[s,n,t] ≥ b_snt[s,n,t]-b⁰_sn[s,n])
    @constraint(model,
        [s in S, n in N, t in [1]],
        b⁻_snt[s,n,t] ≥ b_snt[s,n,t]-b⁰_sn[s,n])

    # Charge and discharge (t>1)
    @constraint(model,
        [s in S, n in N, t in T[T.>1]],
        b⁺_snt[s,n,t] ≥ b_snt[s,n,t]-b_snt[s,n,t-1])
    @constraint(model,
        [s in S, n in N, t in T[T.>1]],
        b⁻_snt[s,n,t] ≥ b_snt[s,n,t]-b_snt[s,n,t-1])

    # Storage capcity
    @constraint(model,
        [s in S, n in N, t in T],
        b_snt[s,n,t]≤b̄_sn[s,n])

    # Storage
    @constraint(model,
        [s in S, n in N],
        b_snt[s,n,1]==b_snt[s,n,T[end]])

    # Ramping limits
    @constraint(model,
        [g in G, n in N, t in T[T.>1]],
        p_gnt[g,n,t]-p_gnt[g,n,t-1]≥r⁺_g[g])
    @constraint(model,
        [g in G, n in N, t in T[T.>1]],
        p_gnt[g,n,t]-p_gnt[g,n,t-1]≤r⁻_g[g])

    # Voltage angles
    @constraint(model,
        [g in G, l in L′, n in N, n′ in N, t in T[T.>1]],
        (θ_nt[n,t] - θ′_nt[n′,t])*B_l[l] == p_gnt[g,n,t]-p_gnt[g,n′,t])

    return model
end

end # module
