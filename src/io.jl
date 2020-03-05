using CSV, JSON, DataFrames

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

"""Loads parameter values for an instance from CSV and JSON files. Reads the following files from `instance_path`.

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
function Params(instance_path::AbstractString)
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

    # Return Params struct
    Params(
        G, G_r, N, L, T, S, κ, C, C̄, τ, τ_t, Q_gn, A_gnt, D_nt, I_g, M_g, C_g,
        r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b0_sn)
end

"""Save object into JSON file.

# Arguments
- `object`
- `output_path::AbstractString`: Full filepath, e.g., `path.json`.
"""
function save_json(object, filepath::AbstractString)
    open(filepath, "w") do io
        JSON.print(io, object)
    end
end
