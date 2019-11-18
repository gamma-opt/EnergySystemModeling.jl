using JuMP # used for mathematical programming
#using Interact # used for enabling the slider
#using Plots
#using Cbc
#using Suppressor
#using Blink
using Dates

function add_var(
        specs::Vector{Bool},
        model_ES::Model,
        Tech::Vector{Int64},
        Nodes::Vector{Int64},
        Lines::Vector{Int64},
        Tsteps::Array{Int64}
        #, Lol::Array{Int64}
    )

    ## Variables applied in any case
    #@variable(model_ES, f[l in Lines, h in Tsteps])                     # transmission decision
    @variable(model_ES, 0 <= g[t in Tech, n in Nodes, h in Tsteps])     # power output of generators
    @variable(model_ES, 0 <= gen_cap[t in Tech, n in Nodes])            # generation capacity (to be fixed)
    #@variable(model_ES, 0 <= trans_cap[l in Lines])                     # transmission capacity (to be fixed)
    @variable(model_ES, 0 <= Lol[n in Nodes, h in Tsteps])              # loss of load introduced   Having this as a common variable, applies for any case
    #@variable(model_ES, 0 <= RES[n in Nodes])                           # Renewables share to be fixed
    #@fix_value(RES, 0.3, force = true)

    # Investment decision
    if specs[1]
        @variable(model_ES, 0 <= gen_inv[t in Tech, n in Nodes])            # generation expansion
        @variable(model_ES, 0 <= trans_inv[l in Lines])                     # transmission expansion
        #@variable(model_ES, 0 <= Gbar[t in Tech])                           # generation expansion budget (to be fixed)
        #@variable(model_ES, 0 <= Tbar)                  # transmission expansion budget (to be fixed)
        #@fix_value(Gbar)
        #@fix_value(Tbar)
    end

    # UC
    if specs[3]
        @variable(model_ES, u[t in Tech, n in Nodes, h in Tsteps], Bin)     # Binary status of generators
        @variable(model_ES, z[l in Lines, h in Tsteps], Bin)                # Binary status of trans. lines
        @variable(model_ES, su[t in Tech, n in Nodes, h in Tsteps], Bin)    # Binary var. for starting up (generation)
        @variable(model_ES, sd[t in Tech, n in Nodes, h in Tsteps], Bin)    # Binary var. for shutting down (generation)
    end

    # Voltage angles
    if specs[5]
        @variable(model_ES, 0 <= θ[n in Nodes, h in Tsteps, d in 1:2])      # Bus angles (from-to)
    end

    if specs[7]
        @variable(model_ES, 0 <= batteryLevel[n in Nodes, h in Tsteps])                       # Battery level variable
        @variable(model_ES, 0 <= batteryInvestment[n in Nodes])             # investment in Battery for each node
        #@variable(model_ES, 0 <= batteryCapacity[n in Nodes])               # Will be fixed
    end

    return model_ES
end

function add_const(
        specs::Vector{Bool},
        model_ES::Model,                        # JuMP model
        Tech::Vector{Int64},                    # Technologies set
        Nodes::Vector{Int64},                   # Nodes set
        Lines::Vector{Int64},                   # Trans. lines set
        Tsteps::Vector{Int64},                  # Time steps set
        g_min::Vector{Float64},                 # Min generation threshold
        B::Vector{Float64},                     # Susceptance
        rup::Vector{Float64},                   # Ramp up limits
        rdw::Vector{Float64},                   # Ramp down limits
        GC::Vector{Float64},                    # Generation cost
        TC::Vector{Float64},                    # Transmission cost
        GEC::Vector{Float64},                   # Generation expansion cost
        TEC::Float64,                           # Transmission expansion cost
        GFC::Vector{Float64},                   # Generation fixed cost
        SUC::Vector{Float64},                   # Starting up cost
        SDC::Vector{Float64},                   # `Shutting`down cost
        TFC::Vector{Float64},                   # Transmission fixed cost
        RTech::Vector{Int64},                   # Renewables listing required for renewables share
        demand::Vector{Float64},                # total demand of each Node at each Tstep
        RES::Float64,                           # renewables share
        Tbar::Float64,                          # transmission max investment
        Gbar::Vector{Float64}                   # generation expansion max investment
        #batteryInvestment:: Array{Float64},    # for investmentlevel of storage
        #batteryLevel::Array{Float64,2},        # battery level at each node and Tsteps
        #f::Vector{Float64}                     # trans flow
    )

    ## Costs
    cost0 = 0
    cost1 = 0
    cost3 = 0
    cost4 = 0

    ## Basic variables
    g = Array{VariableRef,3}(undef,length(Tech),length(Nodes),length(Tsteps))
    f = Array{VariableRef,2}(undef,length(Lines),length(Tsteps))
    gen_cap = Array{VariableRef,2}(undef,length(Tech),length(Nodes))
    trans_cap = Array{VariableRef,1}(undef,length(Lines))
    Lol = Array{VariableRef,2}(undef,length(Nodes),length(Tsteps))
    #RES = Array{VariableRef,1}(undef, 1)

    for h in Tsteps
        for t in Tech
            for n in Nodes
                g[t,n,h] = variable_by_name(model_ES,"g[$t,$n,$h]")
            end
        end
    end

    # for h in Tsteps
    #     for l in Lines
    #         f[l,h] = variable_by_name(model_ES,"f[$l,$h]")
    #     end
    # end

    for t in Tech
        for n in Nodes
            gen_cap[t,n] = variable_by_name(model_ES,"gen_inv[$t,$n]")
        end
    end

    # for l in Lines
    #     trans_cap[l] = variable_by_name(model_ES,"trans_inv[$l]")
    # end

    for h in Tsteps
        for n in Nodes
            Lol[n,h] = variable_by_name(model_ES, "Lol[$n,$h]")
        end
    end

    # for n in Nodes
    #     RES = variable_by_name(model_ES, "RES")
    # end

    cost0 = @expression(model_ES, sum(g[t,n,h]*GC[t] for t in Tech, n in Nodes, h in Tsteps) )#+
    #sum(f[l,h]*TC[l] for l in Lines, h in Tsteps))

    # Investment decision
    if specs[1]
        ## Var definition
        gen_inv = Array{VariableRef,2}(undef,length(Tech),length(Nodes))
        trans_inv = Array{VariableRef,1}(undef,length(Lines))
        #Gbar = Array{VariableRef,1}(undef,length(Tech))
        #Tbar = Array{VariableRef,1}(undef,1)

        #Tbar = variable_by_name(model_ES,"Tbar")
        for t in Tech
            #Gbar[t] = variable_by_name(model_ES,"Gbar[$t]")
            for n in Nodes
                gen_inv[t,n] = variable_by_name(model_ES,"gen_inv[$t,$n]")
            end
        end

        for l in Lines
            trans_inv[l] = variable_by_name(model_ES,"trans_inv[$l]")
        end

        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps], g[t,n,h] <= gen_inv[t,n])                       # Maximum generation
        # @constraint(model_ES, [l in Lines, h in Tsteps], f[l,h] <= trans_inv[l])                         # Maximum transmission
        # @constraint(model_ES, [l in Lines, h in Tsteps], f[l,h] >= -trans_inv[l])                        # Maximum transmission, negative
        @constraint(model_ES, [t in Tech], sum(gen_inv[t,n] for n in Nodes) <= Gbar[t])    # Maximum ge. invest.
        # @constraint(model_ES, sum(trans_inv[l] for l in Lines) <= Tbar)       # Maximum trans. invest.
        # Renewables share    #TODO doesn't work. is now limited by just having a limitation on production
        # @NLconstraint(model_ES, [t in Tech, n in Nodes, h in Tsteps], sum(g[t,n,h] for h in Tsteps, t in RTech)/(sum(g[t,n,h] for h in Tsteps, t in Tech)) >= RES)
        @constraint(model_ES, [t in [1,2,3], n in Nodes, h in Tsteps], sum(g[t,n,h] for h in Tsteps ) <= 600000)

        # Case when battery is introduced
        if specs[7]
            batteryInvestment = Array{VariableRef, 1}(undef, length(Nodes))
            batteryLevel = Array{VariableRef, 2}(undef, length(Nodes), length(Tsteps))
            for n in Nodes
                for h in Tsteps
                    batteryLevel[n,h] = variable_by_name(model_ES,"batteryLevel[$n, $h]")
                end
                batteryInvestment[n] = variable_by_name(model_ES,"batteryInvestment[$n]")
            end

            @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps], batteryLevel[n,h] <= batteryInvestment[n])
            # Storage investments introduced TODO Include both incest and capacity

            @constraint(
                model_ES,
                [t in Tech, n in Nodes, h in Tsteps[Tsteps.>1], l in Lines],
                demand[h] - Lol[n,h] == sum(g[t,n,h] for t in Tech)
                + batteryLevel[n, h-1] - batteryLevel[n, h]
                + sum(Float64[f[l] for l in 1:length(Lines) if Lines[l][2] == n])    #Balance introduced with
                - sum(Float64[f[l] for l in 1:length(Lines) if Lines[l][1] == n])
            )
            #TODO DEMAND

            #NOTE: Outflow for Node 1:
            #sum(Float64[flow[l] for l in 1:length(Lines) if Lines[l][1] == 1])

            #Outflow for Node n:
            #sum(Float64[flow[l] for l in 1:length(Lines) if Lines[l][1] == n])
            #Inflow for Node n:
            #sum(Float64[flow[l] for l in 1:length(Lines) if Lines[l][2] == n])

            batteryInvestmentCost = 40000 #PUT A PROPER NUMBER there TODO do properly

            cost4 = @expression(model_ES, batteryInvestment[n]*batteryInvestmentCost )
        # else
        #     @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps], #, l in Lines],
        #     demand[h] - Lol[n,h] .<= sum(g[t,n,h] for t in Tech)
        #     #+ sum(Float64[f[l] for l in Lines if Lines[l][2] == n])
        #     #- sum(Float64[f[l] for l in Lines if Lines[l][1] == n])
        #     ) #TODO DEMAND
        end

        ## Expansion costs
        cost1 = @expression(model_ES, sum(gen_inv[t,n]*GEC[t] for t in Tech, n in Nodes))# +
        #sum(trans_inv[l]*TEC for l in Lines))
    end

    if specs[3]
        u = Array{VariableRef,3}(undef,length(Tech),length(Nodes),length(Tsteps))
        su = Array{VariableRef,3}(undef,length(Tech),length(Nodes),length(Tsteps))
        sd = Array{VariableRef,3}(undef,length(Tech),length(Nodes),length(Tsteps))
        z = Array{VariableRef,2}(undef,length(Lines),length(Tsteps))

        for h in Tsteps
            for t in Tech
                for n in Nodes
                    u[t,n,h] = variable_by_name(model_ES,"u[$t,$n,$h]")
                    su[t,n,h] = variable_by_name(model_ES,"su[$t,$n,$h]")
                    sd[t,n,h] = variable_by_name(model_ES,"sd[$t,$n,$h]")
                end
            end

            for l in Lines
                z[l,h] = variable_by_name(model_ES,"z[$l,$h]")
            end
        end

        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps], g[t,n,h] <= gen_cap[t,n]*u[t,n,h])           # Generation cap (depending on the investment)
        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps], g[t,n,h] >= g_min[t]*u[t,n,h])               # Minimum start-up power
        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps[Tsteps.>1]], 1 - u[t,n,h-1] >= su[t,n,h])        # If it's commited at t-1 it cannot start up at t
        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps[Tsteps.>1]], u[t,n,h-1] >= sd[t,n,h])            # If it's not commited at t-1 it cannot shut down at t
        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps[Tsteps.>1]], u[t,n,h] - u[t,n,h-1] == su[t,n,h] - sd[t,n,h])     # Plugs u, su, and sd variables

        ## Unit Commitment costs
        cost3 = @expression(model_ES, sum(u[t,n,h].*GFC[t] for t in Tech, n in Nodes, h in Tsteps) +
        sum(su[t,n,h].*SUC[t] for t in Tech, n in Nodes, h in Tsteps) +
        sum(sd[t,n,h].*SDC[t] for t in Tech, n in Nodes, h in Tsteps) +
        sum(z[l,h].*TFC[l] for l in Lines, h in Tsteps) )
    elseif specs[2]
        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps], g[t,n,h] <= gen_cap[t,n])           # Generation cap (depending on the investment)
        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps], g[t,n,h] >= g_min[t])               # Minimum start-up power
    end

    if specs[4]
        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps[Tsteps.>1]], g[t,n,h] - g[t,n,h-1] >= rup[t])
        @constraint(model_ES, [t in Tech, n in Nodes, h in Tsteps[Tsteps.>1]], g[t,n,h-1] - g[t,n,h] >= rdw[t])
    end

    if specs[5]
        θ = Array{VariableRef,3}(undef,length(Nodes),length(Tsteps),2)
        for n in Nodes
            for h in Tsteps
                for i in 1:2
                    θ[n,h,i] = variable_by_name(model_ES,"θ[$n,$h,$i]")
                end
            end
        end

        #Constraint voltage angle
        @constraint(
            model_ES,
            [t in Tech, l in Lines, n in Nodes, n_bar in Nodes, h in Tsteps[Tsteps.>1]],
            ( θ[n,h,1] - θ[n_bar,h,2] ) * B[l] == g[t,n,h] - g[t,n_bar,h]
        )   # Faraday law for accounting voltage angles
    end

    if specs[7]
        @objective(model_ES, Min , cost0 + cost1 + cost3)
    else
        @objective(model_ES, Min , cost0 + cost1 + cost3 + cost4)
    end

    return model_ES
end


# TODO: merge vars and constraints
# TODO: create enumeration class for specs

#### NOTE: specs::Vector{Bool}
# specs[1]: Investment planning (Capacity expansion)
# specs[2]: Economic Dispatch
# specs[3]: Unit Commitment
# specs[4]: Ramping limits / SU & SD costs
# specs[5]: Voltage angles
# specs[6]: Security level
# specs[7]: Storage

specs=[true,false,false,false,false,true]
typeof(specs)

Tech = collect(1:5)
Nodes = collect(1:5)
Lines = collect(1:5)
Tsteps = collect(1:10)
Rtech = [4,5] # INCLUDE HERE RENEWABLE ONES' INDEX
#TODO demand =
demand = 2000*rand(length(Tsteps))#,[2000*rand(length(Tsteps))], [2000*rand(length(Tsteps))], [2000*rand(length(Tsteps))], [2000*rand(length(Tsteps))]]

specs=[true,true,true,true,true,false,false]
g_mini = float.([0,0,0,0,0])
B = float([1,1,1,1,1.1])
ru = float([1,1,1,1,1])
rd = float([1,1,1,1,1])
GenC = float([10,20,30,40,50])
TC = float([5,4,3,2,1])
GEC = float([300,250,200,150,100])
TEC = float(1000)
GFC = float([5,5,5,5,5])
SUC = float([2,2,2,2,2])
SDC = float([0,0,0,0,1])
TFC = float([0,0,0,0,0])
RES = 0.3
Tbar = 500000.0
Gbar = [10000.0,50000.0,100000.0,130000.0,110000.0]
#batteryInvestmentCost = float([])

model_ES = Model()
model_ES = add_var(specs,model_ES,Tech,Nodes,Lines,Tsteps)
model_ES = add_const(specs,model_ES,Tech,Nodes,Lines,Tsteps, g_mini,B,ru,rd,GenC,TC,GEC,TEC,GFC,SUC,SDC,TFC, Rtech, demand, RES, Tbar, Gbar)

# using Gurobi
# optimize!(model_ES, with_optimizer(Gurobi.Optimizer))
