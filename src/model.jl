using Parameters, JuMP

"""Define energy system model type as JuMP.Model."""
const EnergySystemModel = Model

"""Specifation for which constraints to include to the model. Constraints that
are not specified are included by default.

# Arguments
- `renewable_target::Bool`: Whether to include renewables target constraint.
- `storage::Bool`: Whether to include storage constraints.
- `ramping::Bool`: Whether to include ramping constraints.
- `voltage_angles::Bool`: Whether to include voltage angle constraints.
"""
@with_kw struct Specs
    renewable_target::Bool
    carbon_cap::Bool
    nuclear_limit::Bool
    storage::Bool
    ramping::Bool
    voltage_angles::Bool
end

"""Input indices and parameters for the model."""
@with_kw struct Params
    region_n::Array{AbstractString, 1}
    technology_g::Array{AbstractString, 1}
    G::Array{Integer, 1}
    G_r::Array{Integer, 1}
    N::Array{Integer, 1}
    L::Array{Array{Integer, 1}, 1}
    T::Array{Integer, 1}
    S::Array{Integer, 1}
    κ::AbstractFloat
    μ::AbstractFloat 
    C::AbstractFloat
    C̄::AbstractFloat
    C_E::AbstractFloat
    τ::Integer
    τ_t::Array{Integer, 1}
    Q_gn::Array{AbstractFloat, 2}
    Q̄_gn::Array{AbstractFloat, 2}
    A_gnt::Array{AbstractFloat, 3}
    D_nt::Array{AbstractFloat, 2}
    I_g::Array{AbstractFloat, 1}
    M_g::Array{AbstractFloat, 1}
    C_g::Array{AbstractFloat, 1}
    e_g::Array{AbstractFloat, 1}
    E_g::Array{AbstractFloat, 1}
    r⁻_g::Array{AbstractFloat, 1}
    r⁺_g::Array{AbstractFloat, 1}
    I_l::Array{AbstractFloat, 1}
    M_l::Array{AbstractFloat, 1}
    C_l::AbstractFloat
    B_l::AbstractFloat
    e_l::AbstractFloat
    ξ_s::Array{AbstractFloat, 1}
    I_s::Array{AbstractFloat, 1}
    C_s::Array{AbstractFloat, 1}
    b0_sn::Array{AbstractFloat, 2}
    W_nmax::Array{AbstractFloat, 1}
    W_nmin::Array{AbstractFloat, 1}
    f_int::Array{AbstractFloat, 2}
    f′_int::Array{AbstractFloat, 2}
    H_n::Array{AbstractFloat, 1}
    H′_n::Array{AbstractFloat, 1}
    F_onmin::Array{AbstractFloat, 1}
    
end

"""Variable values."""
@with_kw struct Variables
    p_gnt::Array{AbstractFloat, 3}
    p̄_gn::Array{AbstractFloat, 2}
    σ_nt::Array{AbstractFloat, 2}
    f⁺_lt::Array{AbstractFloat, 2}
    f⁻_lt::Array{AbstractFloat, 2}
    f̄_l::Array{AbstractFloat, 1}
    b_snt::Array{AbstractFloat, 3}
    b̄_sn::Array{AbstractFloat, 2}
    b⁺_snt::Array{AbstractFloat, 3}
    b⁻_snt::Array{AbstractFloat, 3}
    θ_nt::Array{AbstractFloat, 2}
    θ′_nt::Array{AbstractFloat, 2}
    w_nt::Array{AbstractFloat, 2}
    f_ont::Array{AbstractFloat, 2}
    f′_ont::Array{AbstractFloat, 2}
    f′′_ont::Array{AbstractFloat, 2}
    h_nt::Array{AbstractFloat, 2}
    h′_nt::Array{AbstractFloat, 2}
end

"""Objective values."""
@with_kw struct Objectives
    f1::AbstractFloat
    f2::AbstractFloat
    f3::AbstractFloat
    f4::AbstractFloat
    f5::AbstractFloat
    f6::AbstractFloat
    f7::AbstractFloat
end

"""Expression values."""
@with_kw struct Expressions
    κ′::AbstractFloat
    μ′::AbstractFloat
    C′_E::AbstractFloat
end

data(a::Number) = a
data(a::JuMP.Containers.DenseAxisArray) = a.data

"""Extract variable values from model.

# Arguments
- `model::EnergySystemModel`
"""
function Variables(model::EnergySystemModel)
    tup = Tuple(value.(model[i]) |> data for i in fieldnames(Variables))
    Variables(tup...)
end

"""Extract objective values from model.

# Arguments
- `model::EnergySystemModel`
"""
function Objectives(model::EnergySystemModel)
    tup = Tuple(value.(model[i]) |> data for i in fieldnames(Objectives))
    Objectives(tup...)
end

"""Compute expression values from the results.

# Arguments
- `parameters::Params`
- `variables::Variables`
"""
function Expressions(parameters::Params, variables::Variables)
    G = parameters.G
    G_r = parameters.G_r
    T = parameters.T
    N = parameters.N
    L = parameters.L

    p_gnt = variables.p_gnt
    h_nt = variables.h_nt
    E_g = parameters.E_g
    e_g = parameters.e_g

    κ′ = (sum(p_gnt[g,n,t] for g in G_r, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T)) /
         (sum(p_gnt[g,n,t] for g in G, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T))
    μ′ = sum(p_gnt[5,n,t] for n in N, t in T) /
        (sum(p_gnt[g,n,t] for g in G, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T))
    C′_E = 1 - ((sum(E_g[g] * sum(p_gnt[g,n,t] for n in N, t in T) / e_g[g] for g in G)) / 1868672938)

    Expressions(κ′, μ′, C′_E)
end

"""Creates the energy system model.

# Arguments
- `parameters::Params`
- `specs::Specs`
"""
function EnergySystemModel(parameters::Params, specs::Specs)
    @unpack region_n, technology_g, G, G_r, N, L, T, S, κ, μ, C, C̄, C_E, τ, τ_t, Q_gn, Q̄_gn, A_gnt, D_nt, I_g, M_g,
            C_g, e_g, E_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, e_l, ξ_s, I_s, C_s, b0_sn,
            W_nmax, W_nmin, f_int, f′_int, H_n, H′_n, F_onmin =
            parameters

    # Indices of lines L
    L′ = 1:length(L)

    # Create an instance of JuMP model.
    model = EnergySystemModel()

    ## -- Variables --
    @variable(model, p_gnt[g in G, n in N, t in T]≥0)
    @variable(model, p̄_gn[g in G, n in N]≥0)
    @variable(model, σ_nt[n in N, t in T]≥0)
    @variable(model, f⁺_lt[l in L′, t in T]≥0)
    @variable(model, f⁻_lt[l in L′, t in T]≥0)
    @variable(model, f̄_l[l in L′])
    @variable(model, b_snt[s in S, n in N, t in T]≥0)
    @variable(model, b̄_sn[s in S, n in N]≥0)
    @variable(model, b⁺_snt[s in S, n in N, t in T]≥0)
    @variable(model, b⁻_snt[s in S, n in N, t in T]≥0)
    @variable(model, θ_nt[n in N, t in T]≥0)
    @variable(model, θ′_nt[n in N, t in T]≥0)
    @variable(model, w_nt[n in N, t in T]≥0)
    @variable(model, f_ont[n in N, t in T]≥0)
    @variable(model, f′_ont[n in N, t in T]≥0)
    @variable(model, f′′_ont[n in N, t in T]≥0)
    @variable(model, h_nt[n in N, t in T]≥0)
    @variable(model, h′_nt[n in N, t in T]≥0)

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
        sum(C_l*(f⁺_lt[l,t] + f⁻_lt[l,t])*τ_t[t] for l in L′, t in T))
    @expression(model, f6,
        sum(I_s[s]*b̄_sn[s,n] for s in S, n in N))
    @expression(model, f7,
        sum(C_s[s]*(b⁺_snt[s,n,t] + b⁻_snt[s,n,t])*τ_t[t] for s in S, n in N, t in T))
    @objective(model, Min, f1 + f2 + f3 + f4 + f5 + f6 + f7)

    ## -- Constraints --
    # Transmission lines to node n
    L⁻(n) = (l for (l,(i,j)) in zip(L′,L) if j==n)
    # Transmission lines from node n
    L⁺(n) = (l for (l,(i,j)) in zip(L′,L) if i==n)

    # Energy balance
    @constraint(model,
        b1[n in N, t in T],
        sum(p_gnt[g,n,t] for g in G) +
        σ_nt[n,t] +
        sum(e_l*f⁺_lt[l,t] - f⁻_lt[l,t] for l in L⁻(n)) -
        sum(f⁺_lt[l,t] - e_l*f⁻_lt[l,t] for l in L⁺(n)) +
        sum(ξ_s[s]*b⁻_snt[s,n,t] - b⁺_snt[s,n,t] for s in S) + 
        h_nt[n,t] ==
        D_nt[n,t])

    # Generation capacity
    @constraint(model,
        g1[g in G, n in N, t in T],
        p_gnt[g,n,t] ≤ A_gnt[g,n,t] * (Q_gn[g,n] + p̄_gn[g,n]))
    @constraint(model,
        g2[g in G, n in N],
        (Q_gn[g,n] + p̄_gn[g,n]) / 1000 ≤ Q̄_gn[g,n] / 1000)

    # Minimum renewables share
    if specs.renewable_target
        @constraint(model, g3,
            ((sum(p_gnt[g,n,t] for g in G_r, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T)) / 1000) ≥
            κ * (sum(p_gnt[g,n,t] for g in G, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T)) / 1000)
    end

    # Maximum nuclear share
    if specs.nuclear_limit
        @constraint(model, g4,
            (sum(p_gnt[5,n,t] for n in N, t in T) / 1000) ≤ (μ * (sum(p_gnt[g,n,t] for g in G, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T)) / 1000))
    end

    #Carbon cap
    if specs.carbon_cap
        @constraint(model, g5,
            (sum(E_g[g] * sum(p_gnt[g,n,t] for n in N, t in T) / e_g[g] for g in G)) / 1000 ≤ (1-C_E) * 1868672938 / 1000)
    end
  
    # Shedding upper bound
    @constraint(model,
        g6[n in N, t in T],
        σ_nt[n,t] ≤ C̄ * D_nt[n,t])

    # Transmission capacity
    @constraint(model,
        t1[l in L′, t in T],
        f⁺_lt[l,t] + f⁻_lt[l,t] ≤ f̄_l[l])

    if specs.storage
        # Storage capacity
        @constraint(model,
            s1[s in S, n in N, t in T],
            b_snt[s,n,t]≤b̄_sn[s,n])

        # Discharge
        @constraint(model,
            s3[s in S, n in N, t in T],
            b⁻_snt[s,n,t]≤b_snt[s,n,t])

        # Charge
        @constraint(model,
            s4[s in S, n in N, t in T],
            ξ_s[s]*b⁺_snt[s,n,t]≤b̄_sn[s,n] - b_snt[s,n,t])

        # Storage levels
        @constraint(model,
            s5[s in S, n in N, t in T[T.>1]],
            b_snt[s,n,t]==b_snt[s,n,t-1] + ξ_s[s]*b⁺_snt[s,n,t-1] - b⁻_snt[s,n,t-1])

        # Storage continuity
        @constraint(model,
            s6[s in S, n in N],
            b_snt[s,n,1]==b_snt[s,n,T[end]])
    end

    if specs.ramping
        # Ramping limits
        @constraint(model,
            r1[g in G, n in N, t in T[T.>1]],
            p_gnt[g,n,t]-p_gnt[g,n,t-1] ≤ r⁺_g[g] * (Q_gn[g,n] + p̄_gn[g,n]))
        @constraint(model,
            r2[g in G, n in N, t in T[T.>1]],
            p_gnt[g,n,t]-p_gnt[g,n,t-1] ≥ -r⁻_g[g] * (Q_gn[g,n] + p̄_gn[g,n]))
    end

    if specs.voltage_angles
        # Voltage angles
        @constraint(model,
            v1[g in G, l in L′, n in N, n′ in N, t in T[T.>1]],
            (θ_nt[n,t] - θ′_nt[n′,t])*B_l[l] == p_gnt[g,n,t]-p_gnt[g,n′,t])
    end

    # Hydro energy
    @constraint(model,
        h1[n in N, t in T],
        W_nmin[n] ≤ w_nt[n,t] ≤ W_nmax[n])
    @constraint(model,
        h2[n in N, t in T[T.>1]],
        w_nt[n,t] == w_nt[n,t-1] + f_int[n,t-1] - f_ont[n,t-1])
    @constraint(model,
        h3[n in N],
        w_nt[n,1] == w_nt[n,T[end]])
    @constraint(model,
        h4[n in N, t in T],
        f_ont[n,t] == f′_ont[n,t] + f′′_ont[n,t])
    @constraint(model,
        h5[n in N, t in T],
        f_ont[n,t] ≥ F_onmin[n])
    @constraint(model,
        h6[n in N, t in T],
        0 ≤ f′_ont[n,t] ≤ H_n[n])
    @constraint(model,
        h7[n in N, t in T],
        0 ≤ h′_nt[n,t] ≤ f′_int[n,t])
    @constraint(model,
        h8[n in N, t in T],
        0 ≤ h′_nt[n,t] ≤ H′_n[n])
    @constraint(model,
        h9[n in N, t in T],
        h_nt[n,t] == h′_nt[n,t] + f′_ont[n,t])
    return model
end
