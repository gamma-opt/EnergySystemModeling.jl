using JuMP
using Gurobi
using EnergySystemModel

# Parameters
specs = Specs(true,true,true,true,true,false,false)
tech = collect(1:5)
nodes = collect(1:5)
lines = collect(1:5)
tsteps = collect(1:10)
rtech = [4,5] # INCLUDE HERE RENEWABLE ONES' INDEX
#TODO demand =
demand = 2000*rand(length(tsteps))#,[2000*rand(length(tsteps))], [2000*rand(length(tsteps))], [2000*rand(length(tsteps))], [2000*rand(length(tsteps))]]

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

model = energy_system_model(
    specs, tech, nodes, lines, tsteps, g_mini, B, ru, rd, GenC, TC, GEC, TEC,
    GFC, SUC, SDC, TFC, rtech, demand, RES, Tbar, Gbar
)

optimize!(model, with_optimizer(Gurobi.Optimizer))
