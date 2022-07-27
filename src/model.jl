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
    transmission::Bool = true
    renewable_target::Bool = false
    carbon_cap::Bool = false
    nuclear_limit::Bool = false
    storage::Bool = false
    ramping::Bool = false
    voltage_angles::Bool = false
    hydro::Bool = false
    hydro_simple::Bool = false
end

"""Input indices and parameters for the model."""
@with_kw struct Params
    region_n::Array{AbstractString, 1}
    max_dem_n::Array{Float64, 1}
    technology_g::Array{AbstractString, 1}
    G::Array{Integer, 1}
    G_r::Array{Integer, 1}
    N::Array{Integer, 1}
    L::Array{Array{Integer, 1}, 1}
    L_ind::Array{Integer, 1}
    T::Array{Integer, 1}
    S::Array{Integer, 1}
    H::Array{Integer, 1}
    κ::AbstractFloat
    μ::AbstractFloat 
    C::AbstractFloat
    C̄::AbstractFloat
    C_E::AbstractFloat
    R_E::AbstractFloat
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
    Wmax_hn::Array{AbstractFloat, 2}
    Wmin_hn::Array{AbstractFloat, 2}
    Hmax_hn::Array{AbstractFloat, 2}
    Hmin_hn::Array{AbstractFloat, 2}
    HRmax_n::Array{AbstractFloat, 1}
    Fmin_n::Array{AbstractFloat, 1}
    AH_nt::Array{AbstractFloat, 2}
    AR_nt::Array{AbstractFloat, 2}
    I_h::Vector{Float64}
    M_h::Vector{Float64}
    C_h::Vector{Float64}
    e_h::Vector{Float64}
    E_h::Vector{Float64}
    r⁻_h::Vector{Float64}
    r⁺_h::Vector{Float64}
end

"""Retrieving data from objects typed JuMP.Containers.DenseAxisArray"""
retrieve_data(a::Number) = a
retrieve_data(a::JuMP.Containers.DenseAxisArray) = (JuMP.value.(a)).data
retrieve_data(a::AffExpr) = JuMP.value.(a)
retrieve_data(a::Vector{AffExpr}) = JuMP.value.(a);

"""Extract objective values from a JuMP model.
# Arguments
- `model::EnergySystemModel`
- `JuMPObjDict::::Union{Dict{String, Float64}, Dict{String, Any}}`
"""
function JuMPObj(model::EnergySystemModel, JuMPObjDict::Dict{String, Any})
    dict = Dict(i => model[Symbol(i)] |> retrieve_data |> first for i in keys(JuMPObjDict))
    return dict
end

"""Extract objective values from a JuMP model.
# Arguments
- `model::EnergySystemModel`
- `JuMPObjDict::::Union{Dict{String, Float64}, Dict{String, Any}}`
"""
function JuMPVar(model::EnergySystemModel, JuMPVarDict::Dict{String, Any})
    dict = Dict(i => model[Symbol(i)] |> retrieve_data for i in keys(JuMPVarDict))
    return dict
end

"""Compute expression values from the results.
# Arguments
- `parameters::Params`
- `variables::Variables`
"""
function Expressions(parameters::Params, specs::Specs, variables::Dict{String, Array{Float64}})
    @unpack G, G_r, N, T, H, R_E, e_g, E_g, τ_t = parameters

    """Expression values:
    # Attributes
    - `κ′: Renewable generation share
    - `μ′: Hydro share
    - `C′_E: CO2 emission reduction
    """
    p_gnt = variables["p_gnt"]

    if (specs.hydro || specs.hydro_simple)
        h_hnt = variables["h_hnt"]
        hr_nt = variables["hr_nt"]

        κ′ = (sum(p_gnt[g,n,t]*τ_t[t] for g in G_r, n in N, t in T) + sum(h_hnt[h,n,t]*τ_t[t] for h in H, n in N, t in T) + sum(hr_nt[n,t]*τ_t[t] for n in N, t in T)) /
            (sum(p_gnt[g,n,t]*τ_t[t] for g in G, n in N, t in T) + sum(h_hnt[h,n,t]*τ_t[t] for h in H, n in N, t in T) + sum(hr_nt[n,t]*τ_t[t] for n in N, t in T))
        μ′ = sum(p_gnt[5,n,t]*τ_t[t] for n in N, t in T) /
            (sum(p_gnt[g,n,t]*τ_t[t] for g in G, n in N, t in T) + sum(h_hnt[h,n,t]*τ_t[t] for h in H, n in N, t in T) + sum(hr_nt[n,t]*τ_t[t] for n in N, t in T))
    else
        κ′ = (sum(p_gnt[g,n,t]*τ_t[t] for g in G_r, n in N, t in T)) /
            (sum(p_gnt[g,n,t]*τ_t[t] for g in G, n in N, t in T))
        μ′ = sum(p_gnt[5,n,t]*τ_t[t] for n in N, t in T) /
            (sum(p_gnt[g,n,t]*τ_t[t] for g in G, n in N, t in T))
    end

    C′_E = 1 - (sum(E_g[g] * (sum(p_gnt[g,n,t]*τ_t[t] for n in N, t in T)) / e_g[g] for g in G)) / R_E

    ExpressionsDict = Dict("κ′" => κ′, "μ′" => μ′, "C′_E" => C′_E)

    return ExpressionsDict
end

"""Creates the energy system model.
# Arguments
- `parameters::Params`
- `specs::Specs`
"""
function EnergySystemModel(parameters::Params, specs::Specs)
    @unpack max_dem_n, G, G_r, N, L, L_ind, T, S, H, κ, μ, C, C̄, C_E, R_E, τ_t, Gmin_gn, Gmax_gn, A_gnt, D_nt, I_g, M_g, 
            C_g, e_g, E_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, e_l, Tmin_l, Tmax_l, ξ_s, I_s, C_s, Smin_sn, Smax_sn,
            Wmax_hn, Wmin_hn, Hmax_hn, Hmin_hn, HRmax_n, Fmin_n, AH_nt, AR_nt,
            I_h, M_h, r⁻_h, r⁺_h =
            parameters

    # TODO: include hydro environment constraints: i.e., use C_h, e_h, E_h

    ## Assumptions and caveats
    # i) Numerical scaling: divisions made through the model (usually per 1000) are done to improve numerical scale

    """Variable values."""
    VariablesDict = Dict{String, Any}()

    """Objective values:
    # Attributes
    - `f1: Generation investment and maintenance
    - `f2: Generation operational cost
    - `f3: Shedding cost
    - `f4: Transmission investment and maintenance
    - `f5: Transmission operational cost
    - `f6: Storage investment cost
    - `f7: Storage operational cost
    - `f8: Hydro investment
    """
    ObjectivesDict = Dict{String, Any}()

    # Indices of lines L: L_ind

    # Create an instance of JuMP model.
    model = EnergySystemModel()

    # -- Main variables --
    @variable(model, p_gnt[g in G, n in N, t in T] ≥ 0)
    VariablesDict["p_gnt"] = p_gnt
    @variable(model, Gmax_gn[g,n] ≥ p̄_gn[g in G, n in N] ≥ Gmin_gn[g,n])
    VariablesDict["p̄_gn"] = p̄_gn
    @variable(model, σ_nt[n in N, t in T] ≥ 0)
    VariablesDict["σ_nt"] = σ_nt

    # Transmission variables
    if specs.transmission
        @variable(model, f_lt[l in L_ind, t in T])
        VariablesDict["f_lt"] = f_lt
        @variable(model, f_abs_lt[l in L_ind, t in T] ≥ 0)
        VariablesDict["f_abs_lt"] = f_abs_lt
        @variable(model, Tmax_l[l] ≥ f̄_l[l in L_ind] ≥ Tmin_l[l])
        VariablesDict["f̄_l"] = f̄_l
    end

    if specs.storage
        # Storage variables
        @variable(model, b_snt[s in S, n in N, t in T] ≥ 0)
        VariablesDict["b_snt"] = b_snt
        @variable(model, Smax_sn[s,n] ≥ b̄_sn[s in S, n in N] ≥ Smin_sn[s,n])
        VariablesDict["b̄_sn"] = b̄_sn
        @variable(model, b⁺_snt[s in S, n in N, t in T] ≥ 0)
        VariablesDict["b⁺_snt"] = b⁺_snt
        @variable(model, b⁻_snt[s in S, n in N, t in T] ≥ 0)
        VariablesDict["b⁻_snt"] = b⁻_snt
    end

    if specs.voltage_angles
        # Voltage variables
        @variable(model, θ_nt[n in N, t in T] ≥ 0)
        VariablesDict["θ_nt"] = θ_nt
        @variable(model, θ′_nt[n in N, t in T] ≥ 0)
        VariablesDict["θ′_nt"] = θ′_nt
    end

    if specs.hydro
        # Hydro energy variables
        @variable(model, w_hnt[h in H, n in N, t in T]  ≥ Wmin_hn[h,n])
        VariablesDict["w_hnt"] = w_hnt
        @variable(model, h_hnt[h in H, n in N, t in T] ≥ 0)
        VariablesDict["h_hnt"] = h_hnt
        @variable(model, hr_nt[n in N, t in T] ≥ 0)
        VariablesDict["hr_nt"] = hr_nt
        @variable(model, h̄_hn[h in H, n in N] ≥ Hmin_hn[h,n])
        VariablesDict["h̄_hn"] = h̄_hn
    end

    if specs.hydro_simple
        if !(specs.hydro)
            # Hydro energy variables
            @variable(model, h_hnt[h in H, n in N, t in T] ≥ 0)
            VariablesDict["h_hnt"] = h_hnt
            @variable(model, h̄_hn[h in H, n in N] ≥ Hmin_hn[h,n])
            VariablesDict["h̄_hn"] = h̄_hn
            @variable(model, 0 ≤ hr_nt[n in N, t in T] ≤ 0)
            VariablesDict["hr_nt"] = hr_nt
        end
        # Compute maximal levels 
        a_n = zeros(length(N))
        for n in N
            a_n[n] = sum(AH_nt[n,:])/(sum(Hmin_hn[:,n]) * length(T))
            isnan(a_n[n]) ? a_n[n] = 0 : a_n[n] = a_n[n]
        end
    end

    ## -- Objective --
    # Investment and maintenance of generation capacity
    @expression(model, f1,
        sum(I_g[g]*(p̄_gn[g,n]-Gmin_gn[g,n]) + M_g[g]*p̄_gn[g,n] for g in G, n in N))
    ObjectivesDict["f1"] = f1
    
    # Operational cost of generation dispatch
    @expression(model, f2,
        sum(C_g[g]*p_gnt[g,n,t]*τ_t[t] for g in G, n in N, t in T))
    ObjectivesDict["f2"] = f2
    
    # Shedding cost
    @expression(model, f3,
        sum(C*σ_nt[n,t]*τ_t[t] for n in N, t in T))
    ObjectivesDict["f3"] = f3

    if specs.transmission
        # Investment and maintenance cost of transmission cpacity
        @expression(model, f4,
            sum(I_l[l]*(f̄_l[l]-Tmin_l[l]) + M_l[l]*f̄_l[l] for l in L_ind))
        ObjectivesDict["f4"] = f4

        # Operational cost of transmission flow
        @expression(model, f5,
            sum(C_l[l]*f_abs_lt[l,t]*τ_t[t] for l in L_ind, t in T))
        ObjectivesDict["f5"] = f5
    end
    
    if specs.storage
        # Investment cost of storage capacity
        @expression(model, f6,
            sum(I_s[s]*(b̄_sn[s,n] - Smin_sn[s,n]) for s in S, n in N))
        ObjectivesDict["f6"] = f6

        # Operational cost of storage
        @expression(model, f7,
            sum(C_s[s]*(b⁺_snt[s,n,t] + b⁻_snt[s,n,t])*τ_t[t] for s in S, n in N, t in T))
        ObjectivesDict["f7"] = f7
    end

    # Investment cost of hydro capacity 
    if (specs.hydro || specs.hydro_simple)
        @expression(model, f8,
            sum(I_h*(h̄_hn[h,n] - Hmin_hn[h,n]) + M_h*h̄_hn[h,n] for h in H, n in N))
        ObjectivesDict["f8"] = f8

        # TODO: add operational costs for hydro generation (i.e., C_h)
    end

    @objective(model, Min, sum(sum(flatten(ObjectivesDict[i])) for i in keys(ObjectivesDict))/10^6)

    ## -- Constraints --
    # Transmission lines to node n
    L⁻(n) = (l for (l,(i,j)) in zip(L_ind,L) if j==n)
    # Transmission lines from node n
    L⁺(n) = (l for (l,(i,j)) in zip(L_ind,L) if i==n)

    # TODO: add hydro generation efficiency dependencies (i.e., e_h)

    # Energy balance (dependent on the features selected)
    if specs.transmission && specs.storage && (specs.hydro || specs.hydro_simple)                   # Trans/Stor/Hydro
        @constraint(model,
        b1[n in N, t in T],
        (sum(p_gnt[g,n,t] for g in G) + σ_nt[n,t] + 
        sum(e_l[l]*f_lt[l,t] for l in L⁻(n)) - sum(e_l[l]*f_lt[l,t] for l in L⁺(n)) +
        sum(ξ_s[s]*b⁻_snt[s,n,t] - b⁺_snt[s,n,t] for s in S) +
        sum(h_hnt[h,n,t] for h in H) + hr_nt[n,t])/1000 
        == max_dem_n[n]*D_nt[n,t]/1000)
    elseif specs.transmission && specs.storage && !(specs.hydro || specs.hydro_simple)            # Trans/Stor
        @constraint(model,
        b1[n in N, t in T],
        (sum(p_gnt[g,n,t] for g in G) + σ_nt[n,t] + 
        sum(e_l[l]*f_lt[l,t] for l in L⁻(n)) - sum(e_l[l]*f_lt[l,t] for l in L⁺(n)) +
        sum(ξ_s[s]*b⁻_snt[s,n,t] - b⁺_snt[s,n,t] for s in S))/1000
        == max_dem_n[n]*D_nt[n,t]/1000)
    elseif specs.transmission && !(specs.storage) && (specs.hydro || specs.hydro_simple)            # Trans/Hydro
        @constraint(model,
        b1[n in N, t in T],
        (sum(p_gnt[g,n,t] for g in G) + σ_nt[n,t] + 
        sum(e_l[l]*f_lt[l,t] for l in L⁻(n)) - sum(e_l[l]*f_lt[l,t] for l in L⁺(n)) +
        sum(h_hnt[h,n,t] for h in H) + hr_nt[n,t] )/1000
        == max_dem_n[n]*D_nt[n,t]/1000)
    elseif specs.transmission && !(specs.storage) && !(specs.hydro || specs.hydro_simple)         # Trans
        @constraint(model,
        b1[n in N, t in T],
        (sum(p_gnt[g,n,t] for g in G) + σ_nt[n,t] + 
        sum(e_l[l]*f_lt[l,t] for l in L⁻(n)) - sum(e_l[l]*f_lt[l,t] for l in L⁺(n)))/1000
        == max_dem_n[n]*D_nt[n,t]/1000)
    elseif !(specs.transmission) && !(specs.storage) && !(specs.hydro || specs.hydro_simple)      # -
        @constraint(model,
        b1[n in N, t in T],
        (sum(p_gnt[g,n,t] for g in G) + σ_nt[n,t])/1000
        == max_dem_n[n]*D_nt[n,t]/1000)
    elseif !(specs.transmission) && specs.storage && !(specs.hydro || specs.hydro_simple)         # Stor
        @constraint(model,
        b1[n in N, t in T],
        (sum(p_gnt[g,n,t] for g in G) + σ_nt[n,t] + 
        sum(ξ_s[s]*b⁻_snt[s,n,t] - b⁺_snt[s,n,t] for s in S))/1000
        == max_dem_n[n]*D_nt[n,t]/1000)
    elseif !(specs.transmission) && !(specs.storage) && (specs.hydro || specs.hydro_simple)         # Hydro
        @constraint(model,
        b1[n in N, t in T],
        (sum(p_gnt[g,n,t] for g in G) + σ_nt[n,t] + 
        sum(h_hnt[h,n,t] for h in H) + hr_nt[n,t] )/1000
        == max_dem_n[n]*D_nt[n,t]/1000)
    elseif !(specs.transmission) && specs.storage && (specs.hydro || specs.hydro_simple)            # Stor/Hydro
        @constraint(model,
        b1[n in N, t in T],
        (sum(p_gnt[g,n,t] for g in G) + σ_nt[n,t] + 
        sum(ξ_s[s]*b⁻_snt[s,n,t] - b⁺_snt[s,n,t] for s in S) +
        sum(h_hnt[h,n,t] for h in H) + hr_nt[n,t])/1000
        == max_dem_n[n]*D_nt[n,t]/1000)
    end

    # Generation capacity
    @constraint(model,
        g1[g in G, n in N, t in T],
        p_gnt[g,n,t] ≤ A_gnt[g,n,t] * p̄_gn[g,n])

    # @constraint(model,
    #     g2[g in G, n in N],
    #     p̄_gn[g,n] ≤ Gmax_gn[g,n])

    # Minimum renewables share
    if specs.renewable_target
        if (specs.hydro || specs.hydro_simple)
            @constraint(model, g3[n in N],
                ((sum(p_gnt[g,n,t]*τ_t[t] for g in G_r, t in T) + sum(h_hnt[h,n,t]*τ_t[t] for h in H, t in T)) + sum(hr_nt[n,t]*τ_t[t] for t in T)) ≥
                κ * (sum(p_gnt[g,n,t]*τ_t[t] for g in G, t in T) + sum(h_hnt[h,n,t]*τ_t[t] for h in H, t in T) + sum(hr_nt[n,t]*τ_t[t] for t in T)))
        else
            @constraint(model, g3[n in N],
                sum(p_gnt[g,n,t]*τ_t[t] for g in G_r, t in T) ≥
                κ * sum(p_gnt[g,n,t]*τ_t[t] for g in G, t in T))
        end
    end

    # Maximum nuclear share
    if specs.nuclear_limit
        #TODO: include a parameter with the index(es) of nuclear sources as a parameter in io.jl
        @constraint(model, g4,
            sum(p_gnt[5,n,t] for n in N, t in T) / 1000 ≤ μ * (sum(p_gnt[g,n,t] for g in G, n in N, t in T) + 
            sum(h_hnt[h,n,t] for h in H, n in N, t in T) + sum(hr_nt[n,t] for n in N, t in T)) / 1000)
    end

    #Carbon cap
    if specs.carbon_cap
        @constraint(model, g5[n in N],
            (sum(E_g[g] * sum(p_gnt[g,n,t]*τ_t[t] for t in T) / e_g[g] for g in G)) / 1000 ≤ (1-C_E) * R_E / 1000)
    end
  
    # Shedding upper bound
    # @constraint(model,
    #     g6[n in N, t in T],
    #     σ_nt[n,t] ≤ C̄ * max_dem_n[n]*D_nt[n,t])

    if specs.transmission
        # Transmission capacity
        @constraint(model,
            t1[l in L_ind, t in T],
            f_lt[l,t]/1000 ≤ f̄_l[l]/1000)

        @constraint(model,
            t2[l in L_ind, t in T],
            f_lt[l,t]/1000 ≥ -f̄_l[l]/1000)

        # Absolute value of f_lt
        @constraint(model,
            t3[l in L_ind, t in T],
            f_abs_lt[l,t]/1000 ≥ f_lt[l,t]/1000)

        @constraint(model,
            t4[l in L_ind, t in T],
            f_abs_lt[l,t]/1000 ≥ -f_lt[l,t]/1000)

        # @constraint(model,
        #     t6[l in L_ind, t in T],
        #     f̄_l[l] ≤ Tmax_l[l])
    end

    if specs.storage
        # Initial storage policy
        ## TODO: pass the battery storage factor via input argument
        @constraint(model,
            s0[s in S, n in N, t in [1]],
            b_snt[s,n,t]/1000 == 0.5*b̄_sn[s,n]/1000)
        # Storage capacity
        @constraint(model,
            s1[s in S, n in N, t in T],
            b_snt[s,n,t]/1000 ≤ b̄_sn[s,n]/1000)
        # Discharge limits (t = 1)
        @constraint(model,
            s2[s in S, n in N, t in [1]],
            τ_t[t]*ξ_s[s]*b⁻_snt[s,n,t]/1000 ≤ b_snt[s,n,t]/1000)
        # Discharge limits (t > 1)
        @constraint(model,
            s3[s in S, n in N, t in T[T.>1]],
            τ_t[t]*ξ_s[s]*b⁻_snt[s,n,t]/1000 ≤ τ_t[t-1]*b_snt[s,n,t-1]/1000)
        # Charge
        @constraint(model,
            s4[s in S, n in N, t in T],
            τ_t[t]*b⁺_snt[s,n,t] ≤ (b̄_sn[s,n] - b_snt[s,n,t]))
        # Storage balance
        @constraint(model,
            s5[s in S, n in N, t in T[T.>1]],
            b_snt[s,n,t] == (b_snt[s,n,t-1] + τ_t[t]*(b⁺_snt[s,n,t] - ξ_s[s]*b⁻_snt[s,n,t])))
        # Storage continuity
        @constraint(model,
            s6[s in S, n in N],
            b_snt[s,n,1] == b_snt[s,n,T[end]])
        # Storage capacity bounds
        # @constraint(model,
        #     s7[s in S, n in N],
        #     b̄_sn[s,n] ≤ Smax_sn[s,n])
    end

    if specs.ramping
        # Ramping limits
        @constraint(model,
            r1[g in G, n in N, t in T[T.>1]],
            (p_gnt[g,n,t]-p_gnt[g,n,t-1]) ≤ r⁺_g[g] * p̄_gn[g,n])
        @constraint(model,
            r2[g in G, n in N, t in T[T.>1]],
            (p_gnt[g,n,t]-p_gnt[g,n,t-1]) ≥ -r⁻_g[g] * p̄_gn[g,n])
    end

    if specs.voltage_angles
        # Voltage angles
        @constraint(model,
            v1[g in G, l in L_ind, n in N, n′ in N, t in T[T.>1]],
            (θ_nt[n,t] - θ′_nt[n′,t])*B_l[l] == p_gnt[g,n,t]-p_gnt[g,n′,t])
    end

    if specs.hydro
        # Hydro energy constraints

        # Full reservoir policy
        @constraint(model, 
            h0[h in H, n in N],
            w_hnt[h,n,1]/1000 == Wmax_hn[h,n]/1000)
        # Maximum reservoir level
        @constraint(model,
            h1[h in H, n in N, t in T],
            w_hnt[h,n,t]/1000 ≤ Wmax_hn[h,n]/1000)
        # Reservoir balance
        @constraint(model,
            h2[n in N, t in T[T.>1]],
            sum(w_hnt[h,n,t] for h in H)/1000 ≤ (sum(w_hnt[h,n,t-1] for h in H) + (AH_nt[n,t-1] - sum(h_hnt[h,n,t-1] for h in H))*τ_t[t])/1000)
        # Reservoir temporal continuity
        @constraint(model,
            h3[n in N, h in H],
            w_hnt[h,n,1]/1000 == w_hnt[h,n,T[end]]/1000)
        # Define the minimum hydro flow possible to be used
        α_nt = zeros(length(N),length(T))
        for n in N, t in T
            AH_nt[n,t] + AR_nt[n,t] ≤ Fmin_n[n] ? α_nt[n,t] = 0 : α_nt[n,t] = AH_nt[n,t] + AR_nt[n,t] - Fmin_n[n]
        end
        # Minimum hydro flow constraint
        @constraint(model,
            h4[n in N, t in T],
            (sum(h_hnt[h,n,t] for h in H) + hr_nt[n,t])/1000 ≤ α_nt[n,t]/1000)
        # Maximum hydro flow (availability)
        @constraint(model,
            h5[n in N, t in T],
            sum(h_hnt[h,n,t] for h in H)/1000 ≤ AH_nt[n,t]/1000)
        # Maximum hydro flow (capacity)
        @constraint(model,
            h6[h in H, n in N, t in T],
            h_hnt[h,n,t]/1000 ≤ h̄_hn[h,n]/1000)
        # Maximum RoR hydro generation
        @constraint(model,
            h7[n in N, t in T],
            hr_nt[n,t]/1000 ≤ AR_nt[n,t]/1000)
        # Maximum hydro generation
        @constraint(model,
            h8[n in N, t in T],
            hr_nt[n,t]/1000 ≤ HRmax_n[n]/1000)
        # Maximum hydro installed capacity
        # @constraint(model,
        #     h9[h in H, n in N],
        #     h̄_hn[h,n] ≤ Hmax_hn[h,n])
        if specs.ramping
            # Hydro ramping up
            @constraint(model,
                hr1[h in H, n in N, t in T[T.>1]],
                (h_hnt[h,n,t]-h_hnt[h,n,t-1])/1000 ≤ (r⁺_h[h] * h̄_hn[h,n])/1000)
            # Hydro ramping down
            @constraint(model,
                hr2[h in H, n in N, t in T[T.>1]],
                (h_hnt[h,n,t]-h_hnt[h,n,t-1])/1000 ≥ (-r⁻_h[h] * h̄_hn[h,n])/1000)            
        end
    end

    if specs.hydro_simple
        if !(specs.hydro)
        # Maximum hydro flow (capacity)
        @constraint(model,
            hs1[h in H, n in N, t in T],
            h_hnt[h,n,t]/1000 ≤ h̄_hn[h,n]/1000)
            if specs.ramping
                # Hydro ramping up
                @constraint(model,
                    hsr1[h in H, n in N, t in T[T.>1]],
                    (h_hnt[h,n,t]-h_hnt[h,n,t-1])/1000 ≤ (r⁺_h[h] * h̄_hn[h,n])/1000)
                # Hydro ramping down
                @constraint(model,
                    hsr2[h in H, n in N, t in T[T.>1]],
                    (h_hnt[h,n,t]-h_hnt[h,n,t-1])/1000 ≥ (-r⁻_h[h] * h̄_hn[h,n])/1000)
            end
        end
        # Maximum hydro flow (accumulated)
        @constraint(model,
        hs2[n in N],
        sum(τ_t[t]*h_hnt[h,n,t] for h in H, t in T)/1000 ≤ a_n[n]*sum(τ_t[t]*h̄_hn[h,n] for h in H, t in T)/1000)
    end

    return model, VariablesDict, ObjectivesDict
end