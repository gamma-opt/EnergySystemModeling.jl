module EnergySystemModel

using JuMP

export Specs, energy_system_model

"""Specs"""
struct Specs
    # TODO: booleans
end

"""Create energy system model."""
function energy_system_model()::Model
    # Create an instance of JuMP model.
    model = Model()

    # TODO: variables
    # TODO: objective
    # TODO: constraints

    return model
end

end # module
