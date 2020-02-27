using JuMP, JSON, CSV, DataFrames

const variables = [:p_gnt, :p̄_gn, :σ_nt, :f_lt, :f_lt_abs, :f̄_l, :b_snt,
                   :b̄_sn, :b⁺_snt, :b⁻_snt, :θ_nt, :θ′_nt]
const objectives = [:f1, :f2, :f3, :f4, :f5, :f6, :f7]

"""Equivalent annual cost (EAC)

# Arguments
* `cost::Real`: Net present value of a project.
* `n::Integer`: Number of payments.
* `r::Real`: Interest rate.
"""
function equivalent_annual_cost(cost::Real, n::Integer, r::Real)
    factor = if iszero(r) n else (1 - (1 + r)^(-n)) / r end
    cost / factor
end

"""Input indices and parameters for the model."""
struct Parameters
    G::Array{Integer}
    G_r::Array{Integer}
    N::Array{Integer}
    L::Array{Array{Integer}}
    T::Array{Integer}
    S::Array{Integer}
    κ::AbstractFloat
    C::AbstractFloat
    C̄::AbstractFloat
    τ::Integer
    τ_t::Array{Integer}
    Q_gn::Array{AbstractFloat, 2}
    A_gnt::Array{AbstractFloat, 3}
    D_nt::Array{AbstractFloat, 2}
    I_g::Array{AbstractFloat}
    M_g::Array{AbstractFloat}
    C_g::Array{AbstractFloat}
    r⁻_g::Array{AbstractFloat}
    r⁺_g::Array{AbstractFloat}
    I_l::Array{AbstractFloat}
    M_l::Array{AbstractFloat}
    C_l::Array{AbstractFloat}
    B_l::Array{AbstractFloat}
    ξ_s::Array{AbstractFloat}
    I_s::Array{AbstractFloat}
    C_s::Array{AbstractFloat}
    b0_sn::Array{AbstractFloat, 2}
end

"""Specifation for which constraints to include to the model. Constraints that
are not specified are included by default.

# Arguments

- `renewable_target::Bool`: Whether to include renewables target constraint.
- `storage::Bool`: Whether to include storage constraints.
- `ramping::Bool`: Whether to include ramping constraints.
- `voltage_angles::Bool`: Whether to include voltage angle constraints.
"""
Base.@kwdef struct Specs
    renewable_target::Bool
    storage::Bool
    ramping::Bool
    voltage_angles::Bool
end

"""Loads parameter values for an instance from CSV and JSON files and returns Parameters type. Reads the following files from `instance_path`.

- `indices.json` with fields `G`, `G_r`, `N`, `L`, `T`, `S`
- `constants.json` with fields `kappa`, `C`, `C_bar`, `r`
- `nodes/` -- Time clustered data from the nodes with fields `Dem_Inc`, `Load_mod`, `Max_Load`, `Avail_Win`, `Avail_Sol`
  - `1.csv`
  - `2.csv`
  - ...
- `technology.csv` with fields `cost`, `lifetime`, `M`, `fuel_cost_1`, `fuel_cost_2`, `r_minus`, `r_plus`
- `transmission.csv` with fields `M`, `cost`, `dist`, `lifetime`, `C`, `B`
- `storage.csv` with fields `xi`, `cost`, `lifetime`, `C`, `b0_1, ..., b0_n`

# Arguments

- `instance_path::AbstractString`: Path to the instance directory.
"""
function load_parameters(instance_path::AbstractString):: Parameters
    # Load indexes and constant parameters
    indices = JSON.parsefile(joinpath(instance_path, "indices.json"))

    # TODO: implement time period clustering: τ, T, τ_t

    # Load indices. Convert JSON values to right types.
    G = indices["G"] |> Array{Int}
    G_r = indices["G_r"] |> Array{Int}
    N = indices["N"] |> Array{Int}
    L = indices["L"] |> Array{Array{Int}}
    τ = 1
    T = 1:indices["T"]
    S = indices["S"] |> Array{Int}

    # Load constant parameters
    constants = JSON.parsefile(joinpath(instance_path, "constants.json"))
    κ = constants["kappa"]
    C = constants["C"]
    C̄ = constants["C_bar"]
    interest_rate = constants["r"]

    # Load time clustered parameters
    τ_t = ones(length(T))
    # TODO: load from file
    Q_gn = zeros(length(G), length(N))
    D_nt = zeros(length(N), length(T))
    A_gnt = ones(length(G), length(N), length(T))
    for n in N
        # Load node values from CSV files.
        df = CSV.read(joinpath(instance_path, "nodes", "$n.csv")) |> DataFrame
        D_nt[n, :] = df.Dem_Inc[1] .* df.Load_mod .* df.Max_Load[1]
        A_gnt[1, n, :] = df.Avail_Win
        A_gnt[2, n, :] = df.Avail_Sol
    end

    # Load technology parameters
    technology = joinpath(instance_path, "technology.csv") |>
        CSV.read |> DataFrame
    I_g = equivalent_annual_cost.(technology.cost, technology.lifetime,
                                  interest_rate)
    M_g = technology.M
    C_g = technology.fuel_cost_1 ./ technology.fuel_cost_2 ./ 1000
    r⁻_g = technology.r_minus
    r⁺_g = technology.r_plus

    # Load transmission parameters
    transmission = joinpath(instance_path, "transmission.csv") |>
        CSV.read |> DataFrame
    M_l = transmission.M
    I_l = equivalent_annual_cost.(transmission.cost .* transmission.dist .+ M_l,
                                  transmission.lifetime, interest_rate)
    C_l = transmission.C
    B_l = transmission.B

    # Load storage parameters
    storage = joinpath(instance_path, "storage.csv") |>
        CSV.read |> DataFrame
    ξ_s = storage.xi
    I_s = equivalent_annual_cost.(storage.cost, storage.lifetime, interest_rate)
    C_s = storage.C
    b0_sn = storage[:, [Symbol("b0_$n") for n in N]] |> Matrix

    # Return Parameters struct
    return Parameters(
        G, G_r, N, L, T, S, κ, C, C̄, τ, τ_t, Q_gn, A_gnt, D_nt, I_g, M_g, C_g,
        r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b0_sn)
end

data(a::Number) = a
data(a::JuMP.Containers.DenseAxisArray) = a.data

"""Save results, including parameter, variable, and objective values. Writes
the following JSON files to `output_path`.

- `specs.json`
- `parameters.json`
- `variables.json`
- `objectives.json`

# Arguments

- `specs`::Specs
- `parameters::Parameters`
- `model::Model`
- `output_path::AbstractString`
"""
function save_results(specs::Specs, parameters::Parameters, model::Model,
                      output_path::AbstractString)
    # Save specs values
    open(joinpath(output_path, "specs.json"), "w") do io
        JSON.print(io, specs)
    end

    # Save parameter values
    open(joinpath(output_path, "parameters.json"), "w") do io
        JSON.print(io, parameters)
    end

    # Save variable values
    d1 = Dict(variable => value.(model[variable]) |> data
              for variable in variables)
    open(joinpath(output_path, "variables.json"), "w") do io
        JSON.print(io, d1)
    end

    # Save objective values
    d2 = Dict(objective => value.(model[objective]) |> data
              for objective in objectives)
    open(joinpath(output_path, "objectives.json"), "w") do io
        JSON.print(io, d2)
    end
end

"""Creates the energy system model. Returns JuMP model.

# Arguments

- `parameters::Parameters`
- `specs::Specs`
"""
function energy_system_model(parameters::Parameters, specs::Specs)::Model
    # Parameters
    G = parameters.G
    G_r = parameters.G_r
    N = parameters.N
    L = parameters.L
    T = parameters.T
    S = parameters.S
    κ = parameters.κ
    C = parameters.C
    C̄ = parameters.C̄
    τ = parameters.τ
    τ_t = parameters.τ_t
    Q_gn = parameters.Q_gn
    A_gnt = parameters.A_gnt
    D_nt = parameters.D_nt
    I_g = parameters.I_g
    M_g = parameters.M_g
    C_g = parameters.C_g
    r⁻_g = parameters.r⁻_g
    r⁺_g = parameters.r⁺_g
    I_l = parameters.I_l
    M_l = parameters.M_l
    C_l = parameters.C_l
    B_l = parameters.B_l
    ξ_s = parameters.ξ_s
    I_s = parameters.I_s
    C_s = parameters.C_s
    b0_sn = parameters.b0_sn

    # Indices of lines L
    L′ = 1:length(L)

    # Create an instance of JuMP model.
    model = Model()

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
        sum(C_g[g]*p_gnt[g,n,t]*τ_t[t] for g in G, n in N, t in T))
    @expression(model, f3,
        sum(C*σ_nt[n,t]*τ_t[t] for n in N, t in T))
    @expression(model, f4,
        sum((I_l[l]+M_l[l])*f̄_l[l] for l in L′))
    @expression(model, f5,
        sum(C_l[l]*f_lt_abs[l,t]*τ_t[t] for l in L′, t in T))
    @expression(model, f6,
        sum(I_s[s]*b̄_sn[s,n] for s in S, n in N))
    @expression(model, f7,
        sum(C_s[s]*(b⁺_snt[s,n,t] + b⁻_snt[s,n,t])*τ_t[t] for s in S, n in N, t in T))
    @objective(model, Min, f1 + f2 + f3 + f4 + f5 + f6 + f7)

    ## -- Constraints --
    # Energy balance (t=1)
    @constraint(model,
        [s in S, n in N, t in [1]],
        sum(p_gnt[g,n,t] for g in G) +
        σ_nt[n,t] +
        sum(f_lt[l,t] for (l,(i,j)) in zip(L′,L) if i==n) -
        sum(f_lt[l,t] for (l,(i,j)) in zip(L′,L) if j==n) +
        ξ_s[s] * b_snt[s,n,t] ==
        D_nt[n,t])

    # Energy balance (t>1)
    @constraint(model,
        [s in S, n in N, t in T[T.>1]],
        sum(p_gnt[g,n,t] for g in G) +
        σ_nt[n,t] +
        sum(f_lt[l,t] for (l,(i,j)) in zip(L′,L) if i==n) -
        sum(f_lt[l,t] for (l,(i,j)) in zip(L′,L) if j==n) +
        ξ_s[s] * (b_snt[s,n,t] - b_snt[s,n,t-1]) ==
        D_nt[n,t])

    # Generation capacity
    @constraint(model,
        [g in G, n in N, t in T],
        p_gnt[g,n,t] ≤ A_gnt[g,n,t] * (Q_gn[g,n] + p̄_gn[g,n]))

    # Minimum renewables share
    if specs.renewable_target
        @constraint(model,
            sum(p_gnt[g,n,t] for g in G_r, n in N, t in T) ≥
            κ * sum(p_gnt[g,n,t] for g in G, n in N, t in T))
    end

    # Shedding upper bound
    @constraint(model,
        [n in N, t in T],
        σ_nt[n,t] ≤ C̄ * D_nt[n,t])

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

    if specs.storage
        # Charge and discharge (t=1)
        @constraint(model,
            [s in S, n in N, t in [1]],
            b⁺_snt[s,n,t] ≥ b_snt[s,n,t]-b0_sn[s,n])
        @constraint(model,
            [s in S, n in N, t in [1]],
            b⁻_snt[s,n,t] ≥ b_snt[s,n,t]-b0_sn[s,n])

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
    end

    if specs.ramping
        # Ramping limits
        @constraint(model,
            [g in G, n in N, t in T[T.>1]],
            p_gnt[g,n,t]-p_gnt[g,n,t-1]≥r⁺_g[g])
        @constraint(model,
            [g in G, n in N, t in T[T.>1]],
            p_gnt[g,n,t]-p_gnt[g,n,t-1]≤r⁻_g[g])
    end

    if specs.voltage_angles
        # Voltage angles
        @constraint(model,
            [g in G, l in L′, n in N, n′ in N, t in T[T.>1]],
            (θ_nt[n,t] - θ′_nt[n′,t])*B_l[l] == p_gnt[g,n,t]-p_gnt[g,n′,t])
    end

    return model
end
