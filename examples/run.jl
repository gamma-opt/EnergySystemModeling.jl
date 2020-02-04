using JuMP
# using Gurobi

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModel

# TODO: parameter values
(G, G_r, N, L, T, S, κ, C, τ, τ_t, Q_nt, A_gnt, D_nt, I_g, M_g,
 C_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b⁰_sn) = load_parameters("instance")

model = energy_system_model(
     G, G_r, N, L, T, S, κ, C, τ, τ_t, Q_nt, A_gnt, D_nt, I_g, M_g,
     C_g, r⁻_g, r⁺_g, I_l, M_l, C_l, B_l, ξ_s, I_s, C_s, b⁰_sn)

# optimize!(model, with_optimizer(Gurobi.Optimizer))
