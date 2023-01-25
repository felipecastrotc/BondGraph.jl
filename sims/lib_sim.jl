include("lib_flow.jl")
include("lib_models.jl")

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Source of effort
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1

# function hqp1i1!(du, u, p, t)
#     Pi = 0.0
#     Po = 0.0
#     Qi, ω, Q1, Q3, P2, P4 = u[:]
#     ϵ, Δs, d, L, C, ρ, μ, γa, γb, di, Δsi, Ii, ca, Ia, H = hqp1i1_unpack(p)
#     τ, τp = p.input.func, p.input.p
#     # Impeller
#     du[1] = (P2 - P4 + (ρ / H) * (γa * Qi - γb * ω) * ω - darcyq_lam(Δsi, di, μ) * Qi / H) * H / Ii
#     # Shaft
#     du[2] = (τ(t, τp...) - ca * ω - ρ * (γa * Qi - γb * ω) * Qi) / Ia
#     # Pipeline
#     # Upstream
#     # Q1
#     du[3] = (Pi - P2 - darcyq_lam(Δs, d, μ) * Q1 / H) * H / L
#     # Downstream
#     # Q3
#     du[4] = (P4 - Po - darcyq_lam(Δs, d, μ) * Q3 / H) * H / L
#     # Pressures
#     du[5] = (Q1 - Qi) / (C * H)           # P2 Upstream inlet
#     du[6] = (Qi - Q3) / (C * H)           # P4 Downstream coupling pump
# end

# function hqp1i1_unpack(p)
#     tmp = ConvConst()
#     return [
#         p.sys.pipe.ϵ,
#         p.sys.pipe.l,
#         p.sys.pipe.d,
#         p.sys.pipe.I,
#         p.sys.pipe.C,
#         p.sys.fluid.ρ,
#         p.sys.fluid.μ,
#         p.sys.pump.γa,
#         p.sys.pump.γb,
#         p.sys.pump.d,
#         p.sys.pump.l,
#         p.sys.pump.I,
#         p.sys.shaft.c,
#         p.sys.shaft.I,
#         tmp.H,
#     ]
# end

function hqp1i1!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    s = p.sys
    Q_p, ω, Q_1, Q_3, P_2, P_4 = u[:]
    # Impeller with laminar friction
    du[1] = impeller_lam(t, [P_2, P_4, Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    du[2] = shaft_1st(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    du[3] = pipe_lam(t, [P_i, P_2, Q_1], s.pipe, s.fluid, s.conv)
    du[5] = pipe_p(t, [Q_1, Q_p], s.pipe, s.fluid, s.conv)
    # Pipeline Downstream
    du[4] = pipe_lam(t, [P_4, P_o, Q_3], s.pipe, s.fluid, s.conv)
    du[6] = pipe_p(t, [Q_p, Q_3], s.pipe, s.fluid, s.conv)
    
end

# -------------------------------------------------------------------
# Model manual -  1 upstream 2 downstream - Source of effort
# (H)ead Q-flow-rate V-valve (P)ipe-2 (I)mpeller-1
function hqvp2i1!(du, u, p, t)
    Pi = 0.0
    Po = 0.0
    Qp, ω, Q1, Q3, P2, P4 = u[:]
    ϵ, Δs, d, L, C, ρ, μ, γa, γb, di, Δsi, Ii, ca, Ia, H = hqp1i1_unpack(p)
    τ, τp = p.input.func, p.input.p
    # Impeller
    du[1] = (P_2 - P_4 + (ρ / H) * (γa * Q_p - γb * ω) * ω - darcyq_lam(Δs_p, d, μ) * Q_p / H - k_sh * ρ * (d / 2)^2 * ω^2 * (1 - Q_p / q_s)^2 / H) * H / Ii
    # Shaft
    du[2] = (τ(t, τp...) - ca * ω - ρ * (γa * Qi - γb * ω) * Qi) / Ia
    # Pipeline
    # Upstream
    # Q1
    du[3] = (Pi - P2 - (darcyq_haal(Q_1, Δs_1, d, ϵ, μ, ρ) + ΔPValve(Q_1, ρ, Cv, N_1)) * (Q_1^2) / H) * H / L
    # Downstream
    # Q3
    du[4] = (P4 - P6 - darcyq_haal(Q_3, Δs_3, d, ϵ, μ, ρ) * (Q_3^2) / H) * H / L
    # Q5
    du[5] = (P6 - Po - darcyq_haal(Q_5, Δs_5, d, ϵ, μ, ρ) * (Q_5^2) / H) * H / L
    # Pressures
    du[6] = (Q1 - Qi) / (C * H)           # P2 Upstream inlet
    du[7] = (Qi - Q3) / (C * H)           # P4 Downstream coupling pump
    du[8] = (Q3 - Q5) / (C * H)           # P6 Downstream coupling pump
end