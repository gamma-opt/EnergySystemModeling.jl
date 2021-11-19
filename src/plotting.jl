using Plots, StatsPlots, LaTeXStrings

techcolors = [:lightblue :cyan :yellow :darkgreen :lime :gray :orange :brown :blue]

function perform_plotting(Plots_specs::Dict{String, Bool}, parameters::Params, variables::Dict{String, Array{Float64}}, 
    objectives::Dict{String, Float64}, expressions::Dict{String, Float64}, plots_output_path::AbstractString)
    
    @info "Plotting"
    ENV["GKSwstype"]="nul"                      # Prevent opening plots windows
    gr()

    if Plots_specs["p1"]
        ## Plotting part 1: Objective function values
        @info "Plotting OF"
        savefig(plot_objective_values(objectives),
                joinpath(plots_output_path, "pdf", "objectives.pdf"))
        savefig(plot_objective_values(objectives),
                joinpath(plots_output_path, "png", "objectives.png"))
    end

    if Plots_specs["p2"]
        ## Plotting part 2: Dispatch and storage
        @info "Plotting dispatch and storage levels"

        for n in parameters.N
            savefig(plot_generation_dispatch(parameters, variables, expressions, n),
                    joinpath(plots_output_path,"pdf","generation_dispatch_n$n.pdf"))
            savefig(plot_generation_capacities(parameters, variables, expressions, n),
                    joinpath(plots_output_path,"pdf","generation_capacities_n$n.pdf"))
            savefig(plot_storage_level(parameters, variables, expressions, n),
                    joinpath(plots_output_path,"pdf","storage_n$n.pdf"))
            savefig(plot_box(parameters, variables, expressions, n),
                joinpath(plots_output_path,"pdf","boxplot$n.pdf"))
            savefig(plot_generation_dispatch(parameters, variables, expressions, n),
                joinpath(plots_output_path,"png","generation_dispatch_n$n.png"))
            savefig(plot_generation_capacities(parameters, variables, expressions, n),
                joinpath(plots_output_path,"png","generation_capacities_n$n.png"))
            savefig(plot_storage_level(parameters, variables, expressions, n),
                joinpath(plots_output_path,"png","storage_n$n.png"))
            savefig(plot_box(parameters, variables, expressions, n),
                joinpath(plots_output_path,"png","boxplot$n.png"))
        end
    end

    if Plots_specs["p3"]
        ## Plotting part 3: Storage capacity
        @info "Plotting storage capacities"

        savefig(plot_storage_capacities(parameters, variables, expressions),
                joinpath(plots_output_path,"pdf","storage_capacities.pdf")
                )
        savefig(plot_storage_capacities(parameters, variables, expressions),
                joinpath(plots_output_path,"png","storage_capacities.png")
                )
    end

    if Plots_specs["p4"]
        ## Plotting part 4: Generation technologies dispatch levels
        @info "Plotting dispatch all"

        savefig(plot_box_all(parameters, variables, expressions),
                joinpath(plots_output_path,"pdf","boxplotall.pdf"))
        savefig(plot_box_all(parameters, variables, expressions),
                joinpath(plots_output_path,"png","boxplotall.png"))
    end

    if Plots_specs["p5"]
        ## Plotting part 5: Generation capacities (including hydro)
        @info "Plotting capacities stacked"

        savefig(plot_generation_capacities_stacked(parameters, variables, expressions),
                    joinpath(plots_output_path,"pdf","generation_capacities_stacked.pdf"))
                savefig(plot_generation_capacities_stacked(parameters, variables, expressions),
                    joinpath(plots_output_path,"png","generation_capacities_stacked.png"))
    end

    if Plots_specs["p6"]
        ## Plotting part 6: Consolidated dispatch vs demand
        @info "Plotting dispatch all (boxplot)"

        savefig(plot_dispatch_bars(parameters, variables, expressions),
                joinpath(plots_output_path,"pdf","dispatchbars.pdf"))
        savefig(plot_dispatch_bars(parameters, variables, expressions),
                joinpath(plots_output_path,"png","dispatchbars.png"))
    end

    if Plots_specs["p7"]
        ## Plotting part 7: Transmission flow (per line)
        @info "Plotting transmission flow"

        for l in 1:length(parameters.L_ind)
            savefig(plot_transmission_flow(parameters, variables, expressions, l),
                    joinpath(plots_output_path,"pdf","transmission_flow_L$l.pdf"))
            savefig(plot_transmission_flow(parameters, variables, expressions, l),
                    joinpath(plots_output_path,"png","transmission_flow_L$l.png"))

        end
    end

    if Plots_specs["p8"]
        ## Plotting part 8: Transmission capacities
        @info "Plotting transmission capacities"

        savefig(plot_transmission_capacities(parameters, variables, expressions),
                joinpath(plots_output_path,"pdf","transmission_capacities.pdf"))
        savefig(plot_transmission_capacities(parameters, variables, expressions),
                joinpath(plots_output_path,"png","transmission_capacities.png"))
    end

    if Plots_specs["p9"]
        ## Plotting part 9: Consolidated transmission flow
        @info "Plotting transmission flow"

        savefig(plot_transmission_bars(parameters, variables, expressions),
                joinpath(plots_output_path,"pdf","transmission_bars.pdf"))
        savefig(plot_transmission_bars(parameters, variables, expressions),
                joinpath(plots_output_path,"png","transmission_bars.png"))
    end

    if Plots_specs["p10"]
        ## Plotting part 10: Lost of load
        @info "Plotting LoL"

        savefig(plot_loss_of_load(parameters, variables, expressions),
                joinpath(plots_output_path,"pdf","loss_of_load.pdf"))
        savefig(plot_loss_of_load(parameters, variables, expressions),
                joinpath(plots_output_path,"png","loss_of_load.png"))
    end

    if Plots_specs["p11"]
        ## Plotting in development
        @info "Plotting stacked dispatch (not ready)"
    end
end

function plot_generation_dispatch(p_gnt, h_hnt, G, n, T, region_n, technology_g, κ, C_E, κ′, C′_E)
    colors = techcolors
    p = Plots.plot(
        legend=:outertopright,
        size=(800, 600),
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
    plot!(p, T, [h_hnt[1, n, t] for t in T],
          color = colors[9],
          alpha=0.3,
          label="hydro")
    return p
end

function plot_balance_stacked(p_gnt, p̄_gn, h_hnt, h̄_hn, HRmax_n, G, n, T, region_n, technology_g, κ, C_E, κ′, C′_E)
    colors = techcolors
    p = Plots.plot(
        legend=:outertopright,
        size=(800, 600),
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
    plot!(p, T, [h_hnt[1, n, t] for t in T],
          color = colors[9],
          alpha=0.3,
          label="hydro")
    plot!(p, T, [h̄_hn[n] + HRmax_n[n] for t in T],
          color = colors[9],
          label="")
    return p
end

function plot_generation_capacities(p̄_gn, h̄_hn, HRmax_n, G, n, region_n, technology_g, κ, C_E, κ′, C′_E)
    StatsPlots.bar([G;9], [[p̄_gn[g, n] for g in G]; [h̄_hn[n] + HRmax_n[n]]],
        xticks=([G;9], [technology_g;"hydro"]),
        ylabel=L"\bar{p}_{g,n}\,\mathrm{[MW]}",
        title = "Generation capacity by technology in $(region_n[n])\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
        titlefontsize = 10,
        color=techcolors,
        alpha=0.7,
        legend=false,
        size = (800,600))
end

function plot_generation_capacities_stacked(p̄_gn, h̄_hn, HRmax_n, N, H, region_n, technology_g, κ, C_E, κ′, C′_E)
    p̄_gn
    H_tot = sum(h̄_hn[h,:] for h in H) .+ HRmax_n
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
        alpha=0.7,
        legend = :outertopright,
        size = (800,600))
end

function plot_transmission_flow(f_lt, f̄_l, l, L, T, region_n, κ, C_E, κ′, C′_E)
    p = Plots.plot(T, [f_lt[l, t] for t in T],
             xlabel=L"t",
             ylabel=L"f_{l,t}\,\mathrm{[MWh]}",
             title = "Hourly transmission flow between $(region_n[L[l][1]]) and $(region_n[L[l][2]])\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
             titlefontsize = 10,
             alpha=0.7,
             legend=false,
             size=(800, 600))
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
    p = StatsPlots.plot(size=(800, 600),
             ylabel=L"\sum_t f_{l,t}\,\mathrm{[MWh]}",
             title = "Transmission flow by line\nRenewables share = $(round(κ′,digits=3)) ≥ $κ\nCO2 reduction = $(round(C′_E,digits=3)) ≥ $C_E",
             titlefontsize = 10,
             xtickfontsize = 5,
             xticks=(L_ind, lines),
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
        xtickfontsize = 5,
        xticks=(L_ind, lines),
        ylabel=L"\bar{f}_l\,\mathrm{[MW]}",
        legend=false,
        alpha = 0.7,
        size=(800, 600),
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
              size=(800, 600))
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
        alpha=0.7,
        size=(800, 600),
        legend = :outertopright)
end

function plot_loss_of_load(σ_nt, N, T, region_n, κ, C_E, κ′, C′_E)
    p = Plots.plot(
        legend=:outertopright,
        size=(800, 600),
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

function plot_box(p_gnt, h_hnt, G, n, region_n, technology_g, κ, C_E, κ′, C′_E)
    colors = techcolors
    p = StatsPlots.plot(
        size=(800, 600),
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
    StatsPlots.boxplot!([9], [h_hnt[1, n,:]],
        color = colors[9],
        alpha=0.7,
        label=false
    )
    return p
end

function plot_box_all(p_gnt, h_hnt, G, H, technology_g, κ, C_E, κ′, C′_E)
    colors = techcolors
    p = StatsPlots.plot(
        size=(800, 600),
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
    h_t = permutedims(sum(sum(h_hnt[h,:,:] for h in H), dims=1))
    StatsPlots.boxplot!([9], [h_t],
        color = colors[9],
        alpha=0.7,
        label=false
    )
    return p
end

function plot_dispatch_bars(max_dem_n, p_gnt, h_hnt, D_nt, N, T, H, region_n, technology_g, κ, C_E, κ′, C′_E)
    p_gn = sum(p_gnt[:,:,t] for t in T)
    h_n = sum(h_hnt[h,:,t] for h in H, t in T)
    D_n = sum(max_dem_n[:].*D_nt[:,t] for t in T)
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
        alpha=0.7,
        size = (800,600),
        legend = :outertopright)
    bar!(N, D_n, bar_width=0.03, fillcolor = :black, label="demand")
end

# Overload functions for signature: parameters::Params, model::Model, ...

"""Plot objective value and individual objective values."""
function plot_objective_values(objectives::Union{Dict{String, Float64}, Dict{String, Any}})
    obj_names = [
    "gen_inv",
    "gen_oc",
    "shed_oc",
    "trans_inv",
    "trans_oc",
    "sto_inv",
    "sto_oc",
    "hyd_inv"
    ]
    ObjNames = Dict(collect(keys(sort(objectives))) .=> obj_names);
    fs = keys(objectives)
    vs = [objectives[i] for i in fs]
    nms = [ObjNames[i] for i in fs]
    title = string("Objective value: ", round(sum(vs), sigdigits = 2))
    bar(vs,
        xlabel="Objective function",
        xticks = ([1:1:8;],nms),
        ylabel="EUR",
        title=title,
        legend=false,
        size = (800,600))
end

"""Plot generation dispatch."""
function plot_generation_dispatch(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}}, n::Integer)
    plot_generation_dispatch(variables["p_gnt"], variables["h_hnt"], parameters.G, n, parameters.T, parameters.region_n, 
                             parameters.technology_g, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end

"""Plot generation capacities."""
function plot_generation_capacities(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}}, n::Integer)
    plot_generation_capacities(variables["p̄_gn"], variables["h̄_hn"], parameters.HRmax_n, parameters.G, n, parameters.region_n,
                               parameters.technology_g, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end


"""Plot stacked generation capacities as a stacked graph."""
function plot_generation_capacities_stacked(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}})
    plot_generation_capacities_stacked(variables["p̄_gn"], variables["h̄_hn"], parameters.HRmax_n, parameters.N, parameters.H,
                       parameters.region_n, parameters.technology_g, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end


"""Plot transmission flow."""
function plot_transmission_flow(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}}, l::Integer)
    plot_transmission_flow(variables["f_lt"], variables["f̄_l"], l, parameters.L, parameters.T, parameters.region_n, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end

"""Plot transmission bars."""
function plot_transmission_bars(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}})
    plot_transmission_bars(variables["f_lt"], parameters.L, parameters.L_ind, parameters.T,
                           parameters.region_n, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end

"""Plot transmission capacities."""
function plot_transmission_capacities(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}})
    plot_transmission_capacities(variables["f̄_l"], parameters.L, parameters.L_ind, parameters.region_n, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end

"""Plot storage level."""
function plot_storage_level(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}}, n::Integer)
    plot_storage_level(variables["b_snt"], variables["b̄_sn"], parameters.S, n, parameters.T, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end

"""Plot storage capacities."""
function plot_storage_capacities(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}})
    plot_storage_capacities(variables["b̄_sn"], parameters.N, parameters.region_n,
                            parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end

"""Plot loss of load."""
function plot_loss_of_load(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}})
    plot_loss_of_load(variables["σ_nt"], parameters.N, parameters.T, parameters.region_n, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end

"""Plot generation box."""
function plot_box(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}}, n::Integer)
    plot_box(variables["p_gnt"], variables["h_hnt"], parameters.G, n, parameters.region_n,
             parameters.technology_g, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end

function plot_box_all(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}})
    plot_box_all(variables["p_gnt"], variables["h_hnt"], parameters.G, parameters.H,
              parameters.technology_g, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end

function plot_dispatch_bars(parameters::Params, variables::Union{Dict{String, Float64}, Dict{String, Array{Float64}}}, expressions::Union{Dict{String, Float64}, Dict{String, Any}})
    plot_dispatch_bars(parameters.max_dem_n, variables["p_gnt"], variables["h_hnt"], parameters.D_nt, parameters.N, parameters.T, parameters.H,
                       parameters.region_n, parameters.technology_g, parameters.κ, parameters.C_E, expressions["κ′"], expressions["C′_E"])
end
