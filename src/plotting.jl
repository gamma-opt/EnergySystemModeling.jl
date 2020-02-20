using JuMP, Plots, LaTeXStrings

function plot_generation_dispatch(p_gnt, p̄_gn, G, n, T)
    p = plot(legend=:outertopright)
    for g in G
        plot!(p, T, [p_gnt[g, n, t] for t in T],
              alpha=0.3,
              xlabel=L"t",
              ylabel=L"p_{g,n,t}\,\mathrm{[MWh]}",
              label="g$g")
        plot!(p, T, [p̄_gn[g, n] for t in T],
              label="")
    end
    return p
end

function plot_generation_capacities(p̄_gn, G, n)
    bar(G, [p̄_gn[g, n] for g in G],
        xlabel=L"g",
        ylabel=L"\bar{p}_{g,n}\,\mathrm{[MW]}",
        legend=false)
end

function plot_transmission_flow(f_lt, f̄_l, l, T)
    p = plot(T, [f_lt[l, t] for t in T],
             xlabel=L"t",
             ylabel=L"f_{l,t}\,\mathrm{[MWh]}",
             alpha=0.3,
             legend=false)
    plot!(p, T, [f̄_l[l] for t in T])
    return p
end

function plot_transmission_capacities(f̄_l, L)
    L′ = 1:length(L)
    bar(L′, [f̄_l[l] for l in L′],
        xlabel=L"l",
        ylabel=L"\bar{f}_l\,\mathrm{[MW]}",
        legend=false)
end

function plot_storage(b_snt, b̄_sn, S, n, T)
    p = plot(legend=false)
    for s in S
        plot!(p, T, [b_snt[s, n, t] for t in T],
              xlabel=L"t",
              ylabel=L"b_{s,n,t}\,\mathrm{[MWh]}",
              alpha=0.3)
        plot!(p, T, [b̄_sn[s, n] for t in T])
    end
    return p
end

function plot_storage_capacities(b̄_sn, S, n)
    bar(S, [b̄_sn[s, n] for s in S],
        xlabel=L"s",
        ylabel=L"b̄_{s,n}\,\mathrm{[MW]}",
        legend=false)
end

function plot_loss_of_load(σ_nt, N, T)
    p = plot(legend=:outertopright)
    for n in N
        plot!(p, T, [σ_nt[n, t] for t in T],
              xlabel="t",
              ylabel=L"σ_{n,t}",
              label="n$n")
    end
    return p
end

# Overload functions for signature: parameters::Parameters, model::Model, ...

"""Plot objective value and individual objective values."""
function plot_objective_values(model::Model)
    fs = ["f$k" for k in 1:7]
    vs = [value.(model[Symbol(f)]) for f in fs]
    title = string("Objective value: ", round(objective_value(model)))
    bar(fs, vs,
        xlabel="Objective function",
        ylabel="EUR",
        title=title,
        legend=false)
end

"""Plot generation dispatch."""
function plot_generation_dispatch(
        parameters::Parameters, model::Model, n::Integer)
    plot_generation_dispatch(
        value.(model[:p_gnt]), value.(model[:p̄_gn]),
        parameters.G, n, parameters.T)
end

"""Plot generation capacities."""
function plot_generation_capacities(
        parameters::Parameters, model::Model, n::Integer)
    plot_generation_capacities(value.(model[:p̄_gn]), parameters.G, n)
end

"""Plot transmission flow."""
function plot_transmission_flow(
        parameters::Parameters, model::Model, l::Integer)
    plot_transmission_flow(
        value.(model[:f_lt]), value.(model[:f̄_l]),
        l, parameters.T)
end

"""Plot transmission capacities."""
function plot_transmission_capacities(parameters::Parameters, model::Model)
    plot_transmission_capacities(value.(model[:f̄_l]), parameters.L)
end

"""Plot storage."""
function plot_storage(parameters::Parameters, model::Model, n::Integer)
    plot_storage(
        value.(model[:b_snt]), value.(model[:b̄_sn]),
        parameters.S, n, parameters.T)
end

"""Plot storage capacities."""
function plot_storage_capacities(
        parameters::Parameters, model::Model, n::Integer)
    plot_storage_capacities(value.(model[:b̄_sn]), parameters.S, n)
end

"""Plot loss of load."""
function plot_loss_of_load(parameters::Parameters, model::Model)
    plot_loss_of_load(value.(model[:σ_nt]), parameters.N, parameters.T)
end
