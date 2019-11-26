module EnergySystemModel

using JuMP

export Specs, energy_system_model

"""Specs"""
struct Specs
    investment_planning::Bool  # Capcity expansion
    economic_dispatch::Bool
    unit_commitment::Bool
    ramping_limits::Bool  # SU & SD costs
    voltage_angles::Bool
    security_level::Bool
    storage::Bool
end

"""Create energy system model."""
function energy_system_model(
        specs::Specs,
        techs::Vector{Int64},   # The set of technologies
        nodes::Vector{Int64},   # The set of nodes
        lines::Vector{Int64},   # The set of Trans. lines
        tsteps::Vector{Int64},  # The set of time steps
        g_min::Vector{Float64}, # Min generation threshold
        B::Vector{Float64},     # Susceptance
        rup::Vector{Float64},   # Ramp up limits
        rdw::Vector{Float64},   # Ramp down limits
        GC::Vector{Float64},    # Generation cost
        TC::Vector{Float64},    # Transmission cost
        GEC::Vector{Float64},   # Generation expansion cost
        TEC::Float64,           # Transmission expansion cost
        GFC::Vector{Float64},   # Generation fixed cost
        SUC::Vector{Float64},   # Starting up cost
        SDC::Vector{Float64},   # Shutting down cost
        TFC::Vector{Float64},   # Transmission fixed cost
        rtech::Vector{Int64},   # Renewables listing required for renewables share
        demand::Vector{Float64},# total demand of each Node at each Tstep
        RES::Float64,           # renewables share
        tbar::Float64,          # transmission max investment
        gbar::Vector{Float64}   # generation expansion max investment
        )::Model
    # Create an instance of JuMP model.
    model = Model()

    ## Variables
    #@variable(model, f[l in lines, h in tsteps])                # transmission decision
    @variable(model, 0 <= g[t in techs, n in nodes, h in tsteps]) # power output of generators
    @variable(model, 0 <= gen_cap[t in techs, n in nodes])        # generation capacity (to be fixed)
    #@variable(model, 0 <= trans_cap[l in lines])                # transmission capacity (to be fixed)
    @variable(model, 0 <= lol[n in nodes, h in tsteps])          # loss of load introduced   Having this as a common variable, applies for any case
    #@variable(model, 0 <= RES[n in nodes])                      # Renewables share to be fixed
    #@fix_value(RES, 0.3, force = true)

    # Investment Planning
    @variable(model, 0 <= gen_inv[t in techs, n in nodes])  # generation expansion
    @variable(model, 0 <= trans_inv[l in lines])           # transmission expansion
    #@variable(model, 0 <= gbar[t in techs])                # generation expansion budget (to be fixed)
    #@variable(model, 0 <= tbar)                           # transmission expansion budget (to be fixed)
    #@fix_value(gbar)
    #@fix_value(tbar)

    # Unit Commitment
    @variable(model, u[t in techs, n in nodes, h in tsteps], Bin)     # Binary status of generators
    @variable(model, z[l in lines, h in tsteps], Bin)                # Binary status of trans. lines
    @variable(model, su[t in techs, n in nodes, h in tsteps], Bin)    # Binary var. for starting up (generation)
    @variable(model, sd[t in techs, n in nodes, h in tsteps], Bin)    # Binary var. for shutting down (generation)

    # Voltage Angles
    @variable(model, 0 <= θ[n in nodes, h in tsteps, d in 1:2])   # Bus angles (from-to)

    # Storage
    @variable(model, 0 <= battery_level[n in nodes, h in tsteps])  # Battery level variable
    @variable(model, 0 <= battery_investment[n in nodes])          # investment in Battery for each node
    #@variable(model, 0 <= battery_capacity[n in nodes])           # Will be fixed

    ## Costs. Required for defining the objective.
    cost0 = 0
    cost1 = 0
    cost3 = 0
    cost4 = 0

    ## Constraints
    cost0 = @expression(model, sum(g[t,n,h]*GC[t] for t in techs, n in nodes, h in tsteps) )#+
    #sum(f[l,h]*TC[l] for l in lines, h in tsteps))

    if specs.investment_planning
        # Maximum generation
        @constraint(model, [t in techs, n in nodes, h in tsteps], g[t,n,h] <= gen_inv[t,n])
        # Maximum transmission
        # @constraint(model, [l in lines, h in tsteps], f[l,h] <= trans_inv[l])
        # Maximum transmission, negative
        # @constraint(model, [l in lines, h in tsteps], f[l,h] >= -trans_inv[l])
        # Maximum ge. invest.
        @constraint(model, [t in techs], sum(gen_inv[t,n] for n in nodes) <= gbar[t])
        # Maximum trans. invest.
        # @constraint(model, sum(trans_inv[l] for l in lines) <= tbar)
        # Renewables share    #TODO doesn't work. is now limited by just having a limitation on production
        # @NLconstraint(model, [t in techs, n in nodes, h in tsteps], sum(g[t,n,h] for h in tsteps, t in rtech)/(sum(g[t,n,h] for h in tsteps, t in techs)) >= RES)
        @constraint(model, [t in [1,2,3], n in nodes, h in tsteps], sum(g[t,n,h] for h in tsteps ) <= 600000)

        ## Expansion costs
        cost1 = @expression(model, sum(gen_inv[t,n]*GEC[t] for t in techs, n in nodes))# +
        #sum(trans_inv[l]*TEC for l in lines))
    end

    # Case when battery is introduced
    if specs.investment_planning && specs.storage
        @constraint(model, [t in techs, n in nodes, h in tsteps], battery_level[n,h] <= battery_investment[n])
        # Storage investments introduced TODO Include both incest and capacity

        @constraint(
            model,
            [t in techs, n in nodes, h in tsteps[tsteps.>1], l in lines],
            demand[h] - lol[n,h] == sum(g[t,n,h] for t in techs)
            + battery_level[n, h-1] - battery_level[n, h]
            + sum(Float64[f[l] for l in 1:length(lines) if lines[l][2] == n])    #Balance introduced with
            - sum(Float64[f[l] for l in 1:length(lines) if lines[l][1] == n])
        )
        #TODO DEMAND

        #NOTE: Outflow for Node 1:
        #sum(Float64[flow[l] for l in 1:length(lines) if lines[l][1] == 1])

        #Outflow for Node n:
        #sum(Float64[flow[l] for l in 1:length(lines) if lines[l][1] == n])
        #Inflow for Node n:
        #sum(Float64[flow[l] for l in 1:length(lines) if lines[l][2] == n])

        # TODO: should be given as a parameter
        battery_investment_cost = 40000 #PUT A PROPER NUMBER there TODO do properly

        cost4 = @expression(model, battery_investment[n]*battery_investment_cost )
    # elseif specs.investment_planning
    #     @constraint(model, [t in techs, n in nodes, h in tsteps], #, l in lines],
    #     demand[h] - lol[n,h] .<= sum(g[t,n,h] for t in techs)
    #     #+ sum(Float64[f[l] for l in lines if lines[l][2] == n])
    #     #- sum(Float64[f[l] for l in lines if lines[l][1] == n])
    #     ) #TODO DEMAND
    end

    if specs.unit_commitment
        # Generation cap (depending on the investment)
        @constraint(model, [t in techs, n in nodes, h in tsteps], g[t,n,h] <= gen_cap[t,n]*u[t,n,h])
        # Minimum start-up power
        @constraint(model, [t in techs, n in nodes, h in tsteps], g[t,n,h] >= g_min[t]*u[t,n,h])
        # If it's commited at t-1 it cannot start up at t
        @constraint(model, [t in techs, n in nodes, h in tsteps[tsteps.>1]], 1 - u[t,n,h-1] >= su[t,n,h])
        # If it's not commited at t-1 it cannot shut down at t
        @constraint(model, [t in techs, n in nodes, h in tsteps[tsteps.>1]], u[t,n,h-1] >= sd[t,n,h])
        # Plugs u, su, and sd variables
        @constraint(model, [t in techs, n in nodes, h in tsteps[tsteps.>1]], u[t,n,h] - u[t,n,h-1] == su[t,n,h] - sd[t,n,h])

        ## Unit Commitment costs
        cost3 = @expression(model, sum(u[t,n,h].*GFC[t] for t in techs, n in nodes, h in tsteps) +
            sum(su[t,n,h].*SUC[t] for t in techs, n in nodes, h in tsteps) +
            sum(sd[t,n,h].*SDC[t] for t in techs, n in nodes, h in tsteps) +
            sum(z[l,h].*TFC[l] for l in lines, h in tsteps) )
    elseif specs.economic_dispatch
        # Generation cap (depending on the investment)
        @constraint(model, [t in techs, n in nodes, h in tsteps], g[t,n,h] <= gen_cap[t,n])
        # Minimum start-up power
        @constraint(model, [t in techs, n in nodes, h in tsteps], g[t,n,h] >= g_min[t])
    end

    if specs.ramping_limits
        @constraint(model, [t in techs, n in nodes, h in tsteps[tsteps.>1]], g[t,n,h] - g[t,n,h-1] >= rup[t])
        @constraint(model, [t in techs, n in nodes, h in tsteps[tsteps.>1]], g[t,n,h-1] - g[t,n,h] >= rdw[t])
    end

    if specs.voltage_angles
        # Faraday law for accounting voltage angles
        @constraint(model,
            [t in techs, l in lines, n in nodes, n_bar in nodes, h in tsteps[tsteps.>1]],
            ( θ[n,h,1] - θ[n_bar,h,2] ) * B[l] == g[t,n,h] - g[t,n_bar,h])
    end

    ## Objective
    @objective(model, Min , cost0 + cost1 + cost3 + cost4)

    return model
end

end # module
