module EnergySystemModeling

include("model.jl")
include("io.jl")
export EnergySystemModel,
    Specs,
    Params,
    perform_Params,
    JuMPObj,
    JuMPVar,
    retrieve_data,
    Expressions,
    equivalent_annual_cost,
    save_json,
    load_json,
    create_nodedata,
    replace_nans,
    getdispatch,
    read_clusters,
    unpack_clusters_features,
    change_time_parameters,
    break_FTR_parameters_in_periods

include("plotting.jl")
export plot_objective_values,
    plot_generation_dispatch,
    plot_generation_capacities,
    plot_generation_capacities_stacked,
    plot_transmission_flow,
    plot_transmission_capacities,
    plot_transmission_bars,
    plot_storage_level,
    plot_storage_capacities,
    plot_loss_of_load,
    plot_box,
    plot_box_all,
    plot_dispatch_bars,
    perform_plotting

include("aggreg.jl")
export SeriesInstance,
    ClustInstance,
    AggregInstance,
    DistUpdate,
    load_series_instance,
    load_clust_instance,
    aggreg1D,
    cdad,
    sorting_custom,
    search_min_dist,
    compute_dist,
    update_marker,
    replace_lines,
    update_clust!,
    update_k!,
    find_clusters!,
    write_clust_instance!,
    read_clust_instance
end # module
