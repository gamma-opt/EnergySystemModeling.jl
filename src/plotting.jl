using Plots, StatsPlots, LaTeXStrings

function plot_generation_dispatch(p_gnt, p̄_gn, h_nt, H_n, H′_n, G, n, T)
    p = plot(
        legend=:outertopright,
        size=(780, 400)
    )
    for g in G
        plot!(p, T, [p_gnt[g, n, t] for t in T],
              alpha=0.3,
              xlabel=L"t",
              ylabel=L"p_{g,n,t}\,\mathrm{[MWh]}",
              label="g$g")
        plot!(p, T, [p̄_gn[g, n] for t in T],
              label="")
    end
    plot!(p, T, [h_nt[n, t] for t in T],
          alpha=0.3,
          label="h")
    plot!(p, T, [H_n[n] + H′_n[n] for t in T],
          label="")
    return p
end

function plot_generation_capacities(p̄_gn, H_n, H′_n, G, n)
    bar([G;9], [[p̄_gn[g, n] for g in G]; [H_n[n] + H′_n[n]]],
        xticks=[G;9],
        xlabel=L"g",
        ylabel=L"\bar{p}_{g,n}\,\mathrm{[MW]}",
        legend=false)
end

function plot_transmission_flow(f_lt, f̄_l, l, T)
    p = plot(T, [f_lt[l, t] for t in T],
             xlabel=L"t",
             ylabel=L"f_{l,t}\,\mathrm{[MWh]}",
             alpha=0.3,
             legend=false,
             size=(780, 400))
    plot!(p, T, [f̄_l[l] for t in T])
    return p
end

function plot_transmission_capacities(f̄_l, L)
    L′ = 1:length(L)
    bar(L′, [f̄_l[l] for l in L′],
        xticks=L′,
        xlabel=L"l",
        ylabel=L"\bar{f}_l\,\mathrm{[MW]}",
        legend=false)
end

function plot_storage_level(b_snt, b̄_sn, S, n, T)
    p = plot(legend=false)
    for s in S
        plot!(p, T, [b_snt[s, n, t] for t in T],
              xlabel=L"t",
              ylabel=L"b_{s,n,t}\,\mathrm{[MWh]}",
              alpha=0.3,
              size=(780, 400))
        plot!(p, T, [b̄_sn[s, n] for t in T])
    end
    return p
end

function plot_storage_capacities(b̄_sn, S, n)
    bar(S, [b̄_sn[s, n] for s in S],
        xticks=S,
        xlabel=L"s",
        ylabel=L"b̄_{s,n}\,\mathrm{[MW]}",
        legend=false)
end

function plot_loss_of_load(σ_nt, N, T)
    p = plot(
        legend=:outertopright,
        size=(780, 400)
    )
    for n in N
        plot!(p, T, [σ_nt[n, t] for t in T],
              xlabel=L"t",
              ylabel=L"\sigma_{n,t}",
              label="n$n")
    end
    return p
end

function plot_box(p_gnt, h_nt, G, n, T)
    p = plot(
        legend=:outertopright,
        size=(780, 400)
    )
    for g in G
        boxplot!([g], [p_gnt[g,n,:]],
              alpha=0.3,
              xlabel=L"t",
              ylabel=L"p_{g,n,t}\,\mathrm{[MWh]}",
              label="g$g")

    end
    boxplot!([9], [h_nt[n,:]],
        alpha=0.3,
        xlabel=L"g",
        ylabel=L"p_{g,n,t}\,\mathrm{[MWh]}",
        label="g9")
    return p
end

# Overload functions for signature: parameters::Params, model::Model, ...

"""Plot objective value and individual objective values."""
function plot_objective_values(objectives::Objectives)
    fs = fieldnames(Objectives)
    vs = [getfield(objectives, field) for field in fs]
    title = string("Objective value: ", round(sum(vs)))
    bar(vs,
        xlabel="Objective function",
        ylabel="EUR",
        title=title,
        legend=false)
end

"""Plot generation dispatch."""
function plot_generation_dispatch(
        parameters::Params, variables::Variables, n::Integer)
    plot_generation_dispatch(
        variables.p_gnt, variables.p̄_gn, variables.h_nt, parameters.H_n, parameters.H′_n, parameters.G, n, parameters.T)
end

"""Plot generation capacities."""
function plot_generation_capacities(
        parameters::Params, variables::Variables, n::Integer)
    plot_generation_capacities(variables.p̄_gn, parameters.H_n, parameters.H′_n, parameters.G, n)
end

"""Plot transmission flow."""
function plot_transmission_flow(
        parameters::Params, variables::Variables, l::Integer)
    plot_transmission_flow(
        variables.f_lt, variables.f̄_l, l, parameters.T)
end

"""Plot transmission capacities."""
function plot_transmission_capacities(
        parameters::Params, variables::Variables)
    plot_transmission_capacities(variables.f̄_l, parameters.L)
end

"""Plot storage level."""
function plot_storage_level(parameters::Params, variables::Variables, n::Integer)
    plot_storage_level(
        variables.b_snt, variables.b̄_sn, parameters.S, n, parameters.T)
end

"""Plot storage capacities."""
function plot_storage_capacities(
        parameters::Params, variables::Variables, n::Integer)
    plot_storage_capacities(variables.b̄_sn, parameters.S, n)
end

"""Plot loss of load."""
function plot_loss_of_load(parameters::Params, variables::Variables)
    plot_loss_of_load(variables.σ_nt, parameters.N, parameters.T)
end

"""Plot generation box."""
function plot_box(
        parameters::Params, variables::Variables, n::Integer)
    plot_box(
        variables.p_gnt, variables.h_nt, parameters.G, n, parameters.T)
end
