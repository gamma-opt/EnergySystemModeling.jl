module EnergySystemModeling

include("model.jl")
include("io.jl")
export EnergySystemModel,
    Specs,
    Params,
    JuMPObj,
    Expressions,
    equivalent_annual_cost,
    save_json,
    load_json,
    create_nodedata,
    replace_nans,
    getdispatch,
    read_clusters,
    unpack_clusters_features

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
    plot_dispatch_bars

include("aggreg.jl")
export SData,
    InputData,
    SeriesInstance,
    ClustInstance,
    AggregInstance,
    DistUpdate,
    load_series_instance,
    load_clust_instance,
    aggreg1D,
    cdad,
    search_min_dist!,
    compute_dist,
    update_marker,
    replace_lines,
    update_clust!,
    update_k!,
    find_clusters!

end # module
