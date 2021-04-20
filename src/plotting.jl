using Plots, StatsPlots, LaTeXStrings

techcolors = [:lightblue :cyan :yellow :darkgreen :lime :gray :orange :brown :blue]

function plot_generation_dispatch(p_gnt, h_nt, G, n, T, region_n, technology_g, κ, C_E, κ′, C′_E)
    colors = techcolors
    p = Plots.plot(
        legend=:outertopright,
        size=(780, 400),
        title = "Hourly generation dispatch by technology in $(region_n[n])\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10
    )
    for g in G
        plot!(p, T, [p_gnt[g, n, t] for t in T],
              color = colors[g],
              alpha=0.3,
              xlabel=L"t",
              ylabel=L"p_{g,n,t}\,\mathrm{[MWh]}",
              label=technology_g[g])
    end
    plot!(p, T, [h_nt[n, t] for t in T],
          color = colors[9],
          alpha=0.3,
          label="hydro")
    return p
end

function plot_balance_stacked(p_gnt, p̄_gn, h_nt, H_n, H′_n, G, n, T, region_n, technology_g, κ, C_E, κ′, C′_E)
    colors = techcolors
    p = Plots.plot(
        legend=:outertopright,
        size=(780, 400),
        title = "Hourly generation dispatch by technology in $(region_n[n])\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10
    )
    for g in G
        plot!(p, T, [p_gnt[g, n, t] for t in T],
              color = colors[g],
              alpha=0.3,
              xlabel=L"t",
              ylabel=L"p_{g,n,t}\,\mathrm{[MWh]}",
              label=technology_g[g])
        plot!(p, T, [p̄_gn[g, n] for t in T],
              color = colors[g],
              label="")
    end
    plot!(p, T, [h_nt[n, t] for t in T],
          color = colors[9],
          alpha=0.3,
          label="hydro")
    plot!(p, T, [H_n[n] + H′_n[n] for t in T],
          color = colors[9],
          label="")
    return p
end

function plot_generation_capacities(p̄_gn, H_n, H′_n, G, n, region_n, technology_g, κ, C_E, κ′, C′_E)
    StatsPlots.bar([G;9], [[p̄_gn[g, n] for g in G]; [H_n[n] + H′_n[n]]],
        xticks=([G;9], [technology_g;"hydro"]),
        ylabel=L"\bar{p}_{g,n}\,\mathrm{[MW]}",
        title = "Generation capacity by technology in $(region_n[n])\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10,
        color=permutedims(techcolors),
        alpha=0.7,
        legend=false)
end
   

function plot_generation_capacities_stacked(p̄_gn, H_n, H′_n, N, region_n, technology_g, κ, C_E, κ′, C′_E)
    p̄_gn
    H_tot = H_n + H′_n
    dispatches = permutedims([p̄_gn; permutedims(H_tot)])
    StatsPlots.groupedbar(dispatches,
        bar_position = :stack,
        bar_width=0.7,
        xticks=(N, region_n),
        ylabel=L"\bar{p}_{g,n}\,\mathrm{[MW]}",
        title = "Generation capacity by region and technology\n Renewables share = $(round(κ′,digits=3)) ≥ $κ\n CO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10,
        color=techcolors,
        labels=permutedims([technology_g;"hydro"]),
        lw = 0,
        alpha=0.7)
end


function plot_transmission_flow(f_lt, f̄_l, l, L, T, region_n, κ, C_E, κ′, C′_E)
    p = Plots.plot(T, [f_lt[l, t] for t in T],
             xlabel=L"t",
             ylabel=L"f_{l,t}\,\mathrm{[MWh]}",
             title = "Hourly transmission flow between $(region_n[L[l][1]]) and $(region_n[L[l][2]])\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
             titlefontsize = 10,
             alpha=0.7,
             legend=false,
             size=(780, 400))
    plot!(p, T, [f̄_l[l] for t in T], color=:green)
    plot!(p, T, [-f̄_l[l] for t in T], color=:green)
    return p
end

function plot_transmission_bars(f_lt, L, L_ind, T, region_n, κ, C_E, κ′, C′_E)
    lines = Array{AbstractString, 1}(undef, length(L))
    for i in L_ind
        lines[i] = "$(region_n[L[i][1]])-$(region_n[L[i][2]])"
    end
    f_l = sum(f_lt[:,t] for t in T)
    p = StatsPlots.plot(size=(780, 400),
             ylabel=L"\sum_t f_{l,t}\,\mathrm{[MWh]}",
             title = "Transmission by line\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
             titlefontsize = 10,
             xticks=(L_ind, lines),
             xrotation = 90,
             legend=false,
             )
    bar!(L_ind, [f_l[l] for l in L_ind], alpha = 0.7, lw = 0)
    return p
end

function plot_transmission_capacities(f̄_l, L, L_ind, region_n, κ, C_E, κ′, C′_E)
    lines = Array{AbstractString, 1}(undef, length(L))
    for i in L_ind
        lines[i] = "$(region_n[L[i][1]])-$(region_n[L[i][2]])"
    end
    StatsPlots.bar(L_ind, [f̄_l[l] for l in L_ind],
        title = "Transmission capacity by line\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10,
        xticks=(L_ind, lines),
        xrotation = 90,
        ylabel=L"\bar{f}_l\,\mathrm{[MW]}",
        legend=false,
        alpha = 0.7,
        lw = 0)
end

function plot_storage_level(b_snt, b̄_sn, S, n, T, κ, C_E, κ′, C′_E)
    p = Plots.plot(legend=:outertopright,
             title = "Hourly storage levels by storage technology\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
             titlefontsize = 10,
        )
    for s in S
        plot!(p, T, [b_snt[s, n, t] for t in T],
              xlabel=L"t",
              ylabel=L"b_{s,n,t}\,\mathrm{[MWh]}",
              labels="battery",
              alpha=0.7,
              size=(780, 400))
        plot!(p, T, [b̄_sn[s, n] for t in T], labels="battery capacity")
    end
    return p
end

function plot_storage_capacities(b̄_sn, N, region_n, κ, C_E, κ′, C′_E)
    # @show b̄_sn, N, region_n
    capacities = permutedims(b̄_sn)
    StatsPlots.groupedbar(capacities,
        bar_position = :stack,
        bar_width=0.7,
        ylabel=L"\bar{p}_{s,n}\,\mathrm{[MWh]}",
        xticks=(N, region_n),
        labels="battery",
        title = "Storage capacity by region\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",  
        titlefontsize = 10,   
        lw = 0,
        alpha=0.7)
end

function plot_loss_of_load(σ_nt, N, T, region_n, κ, C_E, κ′, C′_E)
    p = Plots.plot(
        legend=:outertopright,
        size=(780, 400),
        title = "Loss of load by region\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10
    )
    for n in N
        plot!(p, T, [σ_nt[n, t] for t in T],
              xlabel=L"t",
              ylabel=L"\sigma_{n,t}\,\mathrm{[MWh]}",
              label="$(region_n[n])",
              alpha=0.7)
    end
    return p
end

function plot_box(p_gnt, h_nt, G, n, region_n, technology_g, κ, C_E, κ′, C′_E)
    colors = techcolors
    p = StatsPlots.plot(
        size=(780, 400),
        xticks=([G;9], [technology_g;"hydro"]),
        ylabel=L"p_{g,n,t}\,\mathrm{[MWh]}",
        title = "Hourly generation dispatch by technology in $(region_n[n])\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10
    )
    for g in G
        StatsPlots.boxplot!([g], [p_gnt[g,n,:]],
            color = colors[g],
            alpha=0.7,
            label=false
        )
    end
    StatsPlots.boxplot!([9], [h_nt[n,:]],
        color = colors[9],
        alpha=0.7,
        label=false
    )
    return p
end

function plot_box_all(p_gnt, h_nt, G, technology_g, κ, C_E, κ′, C′_E)
    colors = techcolors
    p = StatsPlots.plot(
        size=(780, 400),
        xticks=([G;9], [technology_g;"hydro"]),
        ylabel=L"\sum_n p_{g,n,t}\,\mathrm{[MWh]}",
        title = "Hourly generation dispatch by technology\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10
    )
    for g in G
        p_gt = permutedims(sum(p_gnt[g, :, :], dims=1))
        StatsPlots.boxplot!([g], [p_gt],
            color = colors[g],
            alpha=0.7,
            label=false
        )
    end
    h_t = permutedims(sum(h_nt[:,:], dims=1))
    StatsPlots.boxplot!([9], [h_t],
        color = colors[9],
        alpha=0.7,
        label=false
    )
    return p
end

function plot_dispatch_bars(p_gnt, h_nt, D_nt, N, T, region_n, technology_g, κ, C_E, κ′, C′_E)
    p_gn = sum(p_gnt[:,:,t] for t in T)
    h_n = sum(h_nt[:,t] for t in T)
    D_n = sum(D_nt[:,t] for t in T)
    dispatches = permutedims([p_gn; permutedims(h_n)])
    groupedbar(dispatches,
        bar_position = :stack,
        bar_width=0.7,
        ylabel=L"\sum_t p_{g,n,t}\,\mathrm{[MWh]}",
        xticks=(N, region_n),
        title = "Generation dispatch by region and technology\n Renewables share = $(round(κ′,digits=3)) ≥ $κ\n CO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10,
        color=techcolors,
        labels=permutedims([technology_g;"hydro"]),
        lw = 0,
        alpha=0.7)
    bar!(N, D_n, bar_width=0.03, fillcolor=repeat([:black], length(N)), label="demand")
end

# Overload functions for signature: parameters::Params, model::Model, ...

"""Plot objective value and individual objective values."""
function plot_objective_values(objectives::Objectives)
    fs = fieldnames(Objectives)
    vs = [getfield(objectives, field) for field in fs]
    title = string("Objective value: ", round(sum(vs), digits = 4))
    bar(vs,
        xlabel="Objective function",
        ylabel="EUR",
        title=title,
        legend=false)
end

"""Plot generation dispatch."""
function plot_generation_dispatch(parameters::Params, variables::Variables, expressions::Expressions, n::Integer)
    plot_generation_dispatch(variables.p_gnt, variables.p̄_gn, variables.h_nt, variables.̄h_n, variables.̄hr_n, parameters.G,
                             n, parameters.T, parameters.region_n, parameters.technology_g, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end

"""Plot generation capacities."""
function plot_generation_capacities(parameters::Params, variables::Variables, expressions::Expressions, n::Integer)
    plot_generation_capacities(variables.p̄_gn, variables.̄h_n, variables.̄hr_n, parameters.G, n, parameters.region_n,
                               parameters.technology_g, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end


"""Plot stacked generation capacities as a stacked graph."""
function plot_generation_capacities_stacked(parameters::Params, variables::Variables, expressions::Expressions)
    plot_generation_capacities_stacked(variables.p̄_gn, variables.̄h_n, variables.̄hr_n, parameters.N,
                       parameters.region_n, parameters.technology_g, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end


"""Plot transmission flow."""
function plot_transmission_flow(parameters::Params, variables::Variables, expressions::Expressions, l::Integer)
    plot_transmission_flow(variables.f_lt, variables.f̄_l, l, parameters.L, parameters.T, parameters.region_n, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end

"""Plot transmission bars."""
function plot_transmission_bars(parameters::Params, variables::Variables, expressions::Expressions)
    plot_transmission_bars(variables.f_lt, parameters.L, parameters.L_ind, parameters.T,
                           parameters.region_n, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end

"""Plot transmission capacities."""
function plot_transmission_capacities(parameters::Params, variables::Variables, expressions::Expressions)
    plot_transmission_capacities(variables.f̄_l, parameters.L, parameters.L_ind, parameters.region_n, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end

"""Plot storage level."""
function plot_storage_level(parameters::Params, variables::Variables, expressions::Expressions, n::Integer)
    plot_storage_level(variables.b_snt, variables.b̄_sn, parameters.S, n, parameters.T, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end

"""Plot storage capacities."""
function plot_storage_capacities(parameters::Params, variables::Variables, expressions::Expressions)
    plot_storage_capacities(variables.b̄_sn, parameters.N, parameters.region_n,
                            parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end

"""Plot loss of load."""
function plot_loss_of_load(parameters::Params, variables::Variables, expressions::Expressions)
    plot_loss_of_load(variables.σ_nt, parameters.N, parameters.T, parameters.region_n, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end

"""Plot generation box."""
function plot_box(parameters::Params, variables::Variables, expressions::Expressions, n::Integer)
    plot_box(variables.p_gnt, variables.h_nt, parameters.G, n, parameters.region_n,
             parameters.technology_g, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end

function plot_box_all(parameters::Params, variables::Variables, expressions::Expressions)
    plot_box_all(variables.p_gnt, variables.h_nt, parameters.G,
              parameters.technology_g, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end

function plot_dispatch_bars(parameters::Params, variables::Variables, expressions::Expressions)
    plot_dispatch_bars(variables.p_gnt, variables.h_nt, parameters.D_nt, parameters.N, parameters.T,
                       parameters.region_n, parameters.technology_g, parameters.κ, parameters.C_E, expressions.κ′, expressions.C′_E)
end
