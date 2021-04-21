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
    hydro::Bool
end

"""Input indices and parameters for the model."""
@with_kw struct Params
    region_n::Array{AbstractString, 1}
    technology_g::Array{AbstractString, 1}
    G::Array{Integer, 1}
    G_r::Array{Integer, 1}
    N::Array{Integer, 1}
    L::Array{Array{Integer, 1}, 1}
    L_ind::Array{Integer, 1}
    T::Array{Integer, 1}
    S::Array{Integer, 1}
    κ::AbstractFloat
    μ::AbstractFloat 
    C::AbstractFloat
    C̄::AbstractFloat
    C_E::AbstractFloat
    R_E::AbstractFloat
    τ::Integer
    τ_t::Array{Integer, 1}
    Gmin_gn::Array{AbstractFloat, 2}
    Gmax_gn::Array{AbstractFloat, 2}
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
    C_l::Array{AbstractFloat, 1}
    B_l::Array{AbstractFloat, 1}
    e_l::Array{AbstractFloat, 1}
    Tmin_l::Array{AbstractFloat, 1}
    Tmax_l::Array{AbstractFloat, 1}
    ξ_s::Array{AbstractFloat, 1}
    I_s::Array{AbstractFloat, 1}
    C_s::Array{AbstractFloat, 1}
    Smin_sn::Array{AbstractFloat, 2}
    Smax_sn::Array{AbstractFloat, 2}
    Wmax_n::Array{AbstractFloat, 1}
    Wmin_n::Array{AbstractFloat, 1}
    Hmax_n::Array{AbstractFloat, 1}
    Hmin_n::Array{AbstractFloat, 1}
    HRcap_n::Array{AbstractFloat, 1}
    Fmin_n::Array{AbstractFloat, 1}
    AH_nt::Array{AbstractFloat, 2}
    AR_nt::Array{AbstractFloat, 2}
    I_h::Array{AbstractFloat, 1}
    M_h::Array{AbstractFloat, 1}
    C_h::Array{AbstractFloat, 1}
    e_h::Array{AbstractFloat, 1}
    E_h::Array{AbstractFloat, 1}
    r⁻_h::Array{AbstractFloat, 1}
    r⁺_h::Array{AbstractFloat, 1}
end

"""Variable values."""
@with_kw struct Variables
    p_gnt::Array{AbstractFloat, 3}
    p̄_gn::Array{AbstractFloat, 2}
    σ_nt::Array{AbstractFloat, 2}
    f_lt::Array{AbstractFloat, 2}
    f̄_l::Array{AbstractFloat, 1}
    b_snt::Array{AbstractFloat, 3}
    b̄_sn::Array{AbstractFloat, 2}
    b⁺_snt::Array{AbstractFloat, 3}
    b⁻_snt::Array{AbstractFloat, 3}
    θ_nt::Array{AbstractFloat, 2}
    θ′_nt::Array{AbstractFloat, 2}
    w_nt::Array{AbstractFloat, 2}
    h_nt::Array{AbstractFloat, 2}
    hr_nt::Array{AbstractFloat, 2}
    h̄_n::Array{AbstractFloat, 1}
end

"""Objective values:
# Attributes
- `f1: Generation investment and maintenance
- `f2: Generation operational cost
- `f3: Shedding cost
- `f4: Transmission investment and maintenance
- `f5: Transmission operational cost
- `f6: Storage investment cost
- `f7: Storage operational cost
"""
@with_kw struct Objectives
    f1::AbstractFloat
    f2::AbstractFloat
    f3::AbstractFloat
    f4::AbstractFloat
    f5::AbstractFloat
    f6::AbstractFloat
    f7::AbstractFloat
    f8::AbstractFloat
end

"""Expression values:
# Attributes
- `κ′: Renewable generation share
- `μ′: Hydro share
- `C′_E: CO2 emission reduction
"""
@with_kw struct Expressions
    κ′::AbstractFloat
    μ′::AbstractFloat
    C′_E::AbstractFloat
end

"""Retrieving data from objects typed JuMP.Containers.DenseAxisArray"""
retrieve_data(a::Number) = a
retrieve_data(a::JuMP.Containers.DenseAxisArray) = a.data

"""Extract variable values from model.
# Arguments
- `model::EnergySystemModel`
"""
function Variables(model::EnergySystemModel)
    tup = Tuple(value.(model[i]) |> retrieve_data for i in fieldnames(Variables))
    Variables(tup...)
end

"""Extract objective values from model.
# Arguments
- `model::EnergySystemModel`
"""
function Objectives(model::EnergySystemModel)
    tup = Tuple(value.(model[i]) |> retrieve_data for i in fieldnames(Objectives))
    Objectives(tup...)
end

"""Compute expression values from the results.
# Arguments
- `parameters::Params`
- `variables::Variables`
"""
function Expressions(parameters::Params, variables::Variables)
    @unpack region_n, technology_g, G, G_r, N, L, L_ind, T, S, κ, μ, C, C̄, C_E, R_E, τ, τ_t, Gmin_gn, Gmax_gn, A_gnt, D_nt, I_g, M_g,
            C_g, e_g, E_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, e_l, Tmin_l, Tmax_l, ξ_s, I_s, C_s, Smin_sn, Smax_sn,
            Wmax_n, Wmin_n, Hmax_n, Hmin_n, HRcap_n, Fmin_n, AH_nt, AR_nt,
            I_h, M_h, C_h, e_h, E_h, r⁻_h, r⁺_h =
            parameters

    p_gnt = variables.p_gnt
    h_nt = variables.h_nt
    E_g = parameters.E_g
    e_g = parameters.e_g

    κ′ = (sum(p_gnt[g,n,t] for g in G_r, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T)) /
         (sum(p_gnt[g,n,t] for g in G, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T))
    μ′ = sum(p_gnt[5,n,t] for n in N, t in T) /
        (sum(p_gnt[g,n,t] for g in G, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T))
    C′_E = 1 - ((sum(E_g[g] * sum(p_gnt[g,n,t] for n in N, t in T) / e_g[g] for g in G)) / R_E)

    Expressions(κ′, μ′, C′_E)
end

"""Creates the energy system model.
# Arguments
- `parameters::Params`
- `specs::Specs`
"""
function EnergySystemModel(parameters::Params, specs::Specs)
    @unpack region_n, technology_g, G, G_r, N, L, L_ind, T, S, κ, μ, C, C̄, C_E, R_E, τ, τ_t, Gmin_gn, Gmax_gn, A_gnt, D_nt, I_g, M_g, 
            C_g, e_g, E_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, e_l, Tmin_l, Tmax_l, ξ_s, I_s, C_s, Smin_sn, Smax_sn,
            Wmax_n, Wmin_n, Hmax_n, Hmin_n, HRcap_n, Fmin_n, AH_nt, AR_nt,
            I_h, M_h, C_h, e_h, E_h, r⁻_h, r⁺_h =
            parameters

    # Indices of lines L: L_ind

    # Create an instance of JuMP model.
    model = EnergySystemModel()

    ## -- Variables --
    @variable(model, p_gnt[g in G, n in N, t in T]≥0)
    @variable(model, p̄_gn[g in G, n in N]≥Gmin_gn[g,n])
    @variable(model, σ_nt[n in N, t in T]≥0)
    @variable(model, f_lt[l in L_ind, t in T])
    @variable(model, f_abs_lt[l in L_ind, t in T]≥0)
    @variable(model, f̄_l[l in L_ind]≥Tmin_l[l])
    @variable(model, b_snt[s in S, n in N, t in T]≥0)
    @variable(model, b̄_sn[s in S, n in N]≥Smin_sn[s,n])
    @variable(model, b⁺_snt[s in S, n in N, t in T]≥0)
    @variable(model, b⁻_snt[s in S, n in N, t in T]≥0)
    @variable(model, θ_nt[n in N, t in T]≥0)
    @variable(model, θ′_nt[n in N, t in T]≥0)
    @variable(model, w_nt[n in N, t in T]≥Wmin_n[n])
    @variable(model, h_nt[n in N, t in T]≥0)
    @variable(model, hr_nt[n in N, t in T]≥0)
    @variable(model, h̄_n[n in N]≥Hmin_n[n])

    ## -- Objective --
    ## Investment and maintenance of generation capacity
    @expression(model, f1,
        sum(I_g[g]*(p̄_gn[g,n]-Gmin_gn[g,n]) + M_g[g]*p̄_gn[g,n] for g in G, n in N))
    
    ## Operational cost of generation dispatch
    @expression(model, f2,
        sum(C_g[g]*p_gnt[g,n,t]*τ_t[t] for g in G, n in N, t in T))
    
    ## Shedding cost
    @expression(model, f3,
        sum(C*σ_nt[n,t]*τ_t[t] for n in N, t in T))

    ##Investment and maintenance cost of transmission cpacity
    @expression(model, f4,
        sum(I_l[l]*(f̄_l[l]-Tmin_l[l]) + M_l[l]*f̄_l[l] for l in L_ind))
    
    ## Operational cost of transmission flow
    @expression(model, f5,
        sum(C_l[l]*f_abs_lt[l,t]*τ_t[t] for l in L_ind, t in T))
    
    ## Investment cost of storage capacity
    @expression(model, f6,
        sum(I_s[s]*(b̄_sn[s,n] - Smin_sn[s,n]) for s in S, n in N))
   
    ## Operational cost of storage
    @expression(model, f7,
        sum(C_s[s]*(b⁺_snt[s,n,t] + b⁻_snt[s,n,t])*τ_t[t] for s in S, n in N, t in T))

    ## Investment cost of hydro capacity 
    # TODO: create an index h for hydro plants
    @expression(model, f8,
        sum(I_h[1]*(h̄_n[n] - Hmin_n[n]) for n in N))
    
    @objective(model, Min, f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8)

    ## -- Constraints --
    # Transmission lines to node n
    L⁻(n) = (l for (l,(i,j)) in zip(L_ind,L) if j==n)
    # Transmission lines from node n
    L⁺(n) = (l for (l,(i,j)) in zip(L_ind,L) if i==n)

    # Energy balance
    @constraint(model,
        b1[n in N, t in T],
        sum(p_gnt[g,n,t] for g in G) +
        σ_nt[n,t] +
        sum(e_l[l]*f_lt[l,t] for l in L⁻(n)) -
        sum(e_l[l]*f_lt[l,t] for l in L⁺(n)) +
        sum(ξ_s[s]*b⁻_snt[s,n,t] - b⁺_snt[s,n,t] for s in S) + 
        h_nt[n,t] +
        hr_nt[n,t] ==
        D_nt[n,t])

    # Generation capacity
    @constraint(model,
        g1[g in G, n in N, t in T],
        p_gnt[g,n,t] ≤ A_gnt[g,n,t] * p̄_gn[g,n])
    @constraint(model,
        g2[g in G, n in N],
        p̄_gn[g,n] ≤ Gmax_gn[g,n])

    # Minimum renewables share
    if specs.renewable_target
        @constraint(model, g3,
            ((sum(p_gnt[g,n,t]*τ_t[t] for g in G_r, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T)) / 1000) ≥
            κ * (sum(p_gnt[g,n,t]*τ_t[t] for g in G, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T)) / 1000)
    end

    # Maximum nuclear share
    if specs.nuclear_limit
        @constraint(model, g4,
            (sum(p_gnt[5,n,t] for n in N, t in T) / 1000) ≤ (μ * (sum(p_gnt[g,n,t] for g in G, n in N, t in T) + sum(h_nt[n,t] for n in N, t in T)) / 1000))
    end

    #Carbon cap
    if specs.carbon_cap
        @constraint(model, g5,
            (sum(E_g[g] * sum(p_gnt[g,n,t] for n in N, t in T) / e_g[g] for g in G)) / 1000 ≤ (1-C_E) * R_E / 1000)
    end
  
    # Shedding upper bound
    @constraint(model,
        g6[n in N, t in T],
        σ_nt[n,t] ≤ C̄ * D_nt[n,t])

    # Transmission capacity
    @constraint(model,
        t1[l in L_ind, t in T],
        f_lt[l,t] ≤ f̄_l[l])

    @constraint(model,
        t2[l in L_ind, t in T],
        f_lt[l,t] ≥ -f̄_l[l])

    # Absolute value of f_lt
    @constraint(model,
        t3[l in L_ind, t in T],
        f_abs_lt[l,t] ≥ f_lt[l,t])

    @constraint(model,
        t4[l in L_ind, t in T],
        f_abs_lt[l,t] ≥ -f_lt[l,t])

    @constraint(model,
        t6[l in L_ind, t in T],
        f̄_l[l] ≤ Tmax_l[l])

    if specs.storage
    ## TODO: revise storage

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

        # Storage capacity bounds
        @constraint(model,
            s7[s in S, n in N],
            b̄_sn[s,n] ≤ Smax_sn[s,n])
    end

    if specs.ramping
        # Ramping limits
        @constraint(model,
            r1[g in G, n in N, t in T[T.>1]],
            p_gnt[g,n,t]-p_gnt[g,n,t-1] ≤ r⁺_g[g] * p̄_gn[g,n])
        @constraint(model,
            r2[g in G, n in N, t in T[T.>1]],
            p_gnt[g,n,t]-p_gnt[g,n,t-1] ≥ -r⁻_g[g] * p̄_gn[g,n])
        if specs.hydro
            @constraint(model,
                r3[g in G, n in N, t in T[T.>1]],
                h_nt[n,t]-h_nt[n,t-1] ≤ r⁺_h[1] * h̄_n[n])

            @constraint(model,
                r4[n in N, t in T[T.>1]],
                h_nt[n,t]-h_nt[n,t-1] ≥ -r⁻_h[1] * h̄_n[n])            
        end
    end

    if specs.voltage_angles
        # Voltage angles
        @constraint(model,
            v1[g in G, l in L_ind, n in N, n′ in N, t in T[T.>1]],
            (θ_nt[n,t] - θ′_nt[n′,t])*B_l[l] == p_gnt[g,n,t]-p_gnt[g,n′,t])
    end

    if specs.hydro
        # Hydro energy
        @constraint(model,
            h1[n in N, t in T],
            w_nt[n,t] ≤ Wmax_n[n])
        @constraint(model,
            h2[n in N, t in T[T.>1]],
            w_nt[n,t] == w_nt[n,t-1] + AH_nt[n,t-1] - h_nt[n,t-1]*τ_t[t])
        @constraint(model,
            h3[n in N],
            w_nt[n,1] == w_nt[n,T[end]])
        @constraint(model,
            h4[n in N, t in T],
            h_nt[n,t] ≤ h̄_n[n] - Fmin_n[n]*AH_nt[n,t])
        @constraint(model,
            h5[n in N, t in T],
            h_nt[n,t] ≤ AH_nt[n,t])
        @constraint(model,
            h6[n in N, t in T],
            h̄_n[n] ≤ Hmax_n[n])
        @constraint(model,
            h7[n in N, t in T],
            hr_nt[n,t] ≤ AR_nt[n,t])
        @constraint(model,
            h8[n in N, t in T],
            hr_nt[n,t] ≤ HRcap_n[n])
    end

    return model
end