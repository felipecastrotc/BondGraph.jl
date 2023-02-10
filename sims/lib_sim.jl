include("lib_flow.jl")
include("lib_models.jl")

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Source of effort
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1
# Legacy version

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

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Source of effort
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1

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
# Model manual -  1 upstream 1 downstream - Source of effort
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1 - Turbulent - Haaland
function hqp1i1th!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    s = p.sys
    Q_p, ω, Q_1, Q_3, P_2, P_4 = u[:]
    # Impeller with laminar friction
    du[1] = impeller_lam(t, [P_2, P_4, Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    du[2] = shaft_1st(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    du[3] = pipe_haal(t, [P_i, P_2, Q_1], s.pipe, s.fluid, s.conv)
    # du[3] = pipe_cheng(t, [P_i, P_2, Q_1], s.pipe, s.fluid, s.conv)
    du[5] = pipe_p(t, [Q_1, Q_p], s.pipe, s.fluid, s.conv)
    # Pipeline Downstream
    du[4] = pipe_haal(t, [P_4, P_o, Q_3], s.pipe, s.fluid, s.conv)
    # du[4] = pipe_cheng(t, [P_4, P_o, Q_3], s.pipe, s.fluid, s.conv)
    du[6] = pipe_p(t, [Q_p, Q_3], s.pipe, s.fluid, s.conv)
end

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Source of effort
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1 - Turbulent - Cheng
function hqp1i1tc!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    s = p.sys
    Q_p, ω, Q_1, Q_3, P_2, P_4 = u[:]
    # Impeller with laminar friction
    du[1] = impeller_lam(t, [P_2, P_4, Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    du[2] = shaft_1st(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    du[3] = pipe_cheng(t, [P_i, P_2, Q_1], s.pipe, s.fluid, s.conv)
    du[5] = pipe_p(t, [Q_1, Q_p], s.pipe, s.fluid, s.conv)
    # Pipeline Downstream
    du[4] = pipe_cheng(t, [P_4, P_o, Q_3], s.pipe, s.fluid, s.conv)
    du[6] = pipe_p(t, [Q_p, Q_3], s.pipe, s.fluid, s.conv)
end

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Source of effort
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1 - Turbulent - Cheng - Shock
function hqp1i1tcs!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    s = p.sys
    Q_p, ω, Q_1, Q_3, P_2, P_4 = u[:]
    # Impeller with laminar friction
    du[1] = impeller_lshock(t, [P_2, P_4, Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    du[2] = shaft_1st(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    du[3] = pipe_cheng(t, [P_i, P_2, Q_1], s.pipe, s.fluid, s.conv)
    du[5] = pipe_p(t, [Q_1, Q_p], s.pipe, s.fluid, s.conv)
    # Pipeline Downstream
    du[4] = pipe_cheng(t, [P_4, P_o, Q_3], s.pipe, s.fluid, s.conv)
    du[6] = pipe_p(t, [Q_p, Q_3], s.pipe, s.fluid, s.conv)
end

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Source of effort
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1 - Turbulent - Cheng - Shock - Valve
function hqvp1i1tcs!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    a = 0.1
    s = p.sys
    Q_p, ω, Q_1, Q_3, P_2, P_4 = u[:]
    # Impeller with laminar friction
    du[1] = impeller_lshock(t, [P_2, P_4, Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    du[2] = shaft_1st(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    du[3] = pipe_cheng(t, [P_i, P_2, Q_1], s.pipe, s.fluid, s.conv)
    du[5] = pipe_p(t, [Q_1, Q_p], s.pipe, s.fluid, s.conv)
    # Pipeline Downstream
    du[4] = pipe_cheng_valve(t, [P_4, P_o, Q_3, a], s.pipe, s.fluid, s.conv)
    du[6] = pipe_p(t, [Q_p, Q_3], s.pipe, s.fluid, s.conv)
end

# -------------------------------------------------------------------
# Model manual -  1 upstream 2 downstream - Source of effort
# (H)ead Q-flow-rate V-valve (P)ipe-2 (I)mpeller-1 - turbulent
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

# -------------------------------------------------------------------
# Model manual -  1 upstream 2 downstream - Source of effort
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1 - Turbulent - Cheng - Shock - Valve
function hqvp2i1tcs!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    a = 0.9
    s = p.sys
    Q_p, ω, Q_1, Q_3, Q_5, P_2, P_4, P_6 = u[:]
    # Impeller with laminar friction
    du[1] = impeller_lshock(t, [P_2, P_4, Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    du[2] = shaft_1st(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    du[3] = pipe_cheng(t, [P_i, P_2, Q_1], s.pipe[1], s.fluid, s.conv)
    du[6] = pipe_p(t, [Q_1, Q_p], s.pipe[1], s.fluid, s.conv)
    # Pipeline downstream
    du[4] = pipe_cheng(t, [P_4, P_6, Q_3], s.pipe[2], s.fluid, s.conv)
    du[7] = pipe_p(t, [Q_p, Q_3], s.pipe[2], s.fluid, s.conv)
    # Pipeline valve downstream
    du[5] = pipe_cheng_valve(t, [P_6, P_o, Q_5, a], s.pipe[3], s.fluid, s.conv)
    du[8] = pipe_p(t, [Q_3, Q_5], s.pipe[3], s.fluid, s.conv)
end


# -------------------------------------------------------------------
# Model manual -  1 upstream 2 downstream - Twin pump
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1 - Turbulent - Cheng - Shock - Valve 
function hqvp2i1tcsc!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    a = 1.0
    s = p.sys
    Q_p, ω, Q_1, Q_3, Q_5, P_2, P_4, P_6 = u[:]
    # Impeller with laminar friction
    # du[1] = impeller_fit1(t, [P_2, P_4, Q_p, ω], s.pump, s.fluid, s.conv)
    du[1] = impeller_lshock(t, [P_2, P_4, Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    du[2] = shaft_1st(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    du[3] = pipe_cheng_twin(t, [P_i, P_2, Q_1], s.pipe[1], s.fluid, s.conv)
    # du[3] = pipe_cheng_cavity(t, [P_i, P_2, Q_1], s.pipe[1], s.fluid, s.conv)
    du[6] = pipe_p(t, [Q_1, Q_p], s.pipe[1], s.fluid, s.conv)
    # Pipeline downstream
    du[4] = pipe_cheng(t, [P_4, P_6, Q_3], s.pipe[2], s.fluid, s.conv)
    du[7] = pipe_p(t, [Q_p, Q_3], s.pipe[2], s.fluid, s.conv)
    # Pipeline valve downstream
    du[5] = pipe_cheng_valve_fit(t, [P_6, P_o, Q_5, a], s.pipe[3], s.fluid, s.conv)
    du[8] = pipe_p(t, [Q_3, Q_5], s.pipe[3], s.fluid, s.conv)
end

# -------------------------------------------------------------------
# Model manual -  3 upstream 3 downstream - Twin pump
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1 - Turbulent - Cheng - Shock - Valve 

function hqvp3i1tcvfit!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    s = p.sys
    Q_p, ω = view(u, 1:2)
    dQ, Q = view(du, 3:8), view(u, 3:8)
    dP, P = view(du, 9:14), view(u, 9:14)
    # Impeller with laminar friction
    du[1] = impeller_fit1(t, [P[3], P[4], Q_p, ω], s.pump, s.fluid, s.conv)
    # du[1] = impeller_lshock(t, [P[3], P[4], Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    # du[2] = shaft_1st(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    du[2] = shaft_1st_fit(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    dQ[1] = pipe_cheng_twin_fit(t, [P_i, P[1], Q[1]], s.pipe[1], s.fluid, s.conv)
    dQ[2] = pipe_cheng(t, [P[1], P[2], Q[2]], s.pipe[2], s.fluid, s.conv)
    dQ[3] = pipe_cheng(t, [P[2], P[3], Q[3]], s.pipe[3], s.fluid, s.conv)
    # Pipeline downstream
    dQ[4] = pipe_cheng(t, [P[4], P[5], Q[4]], s.pipe[4], s.fluid, s.conv)
    # dQ[5] = pipe_lam_valve(t, [P[5], P[6], Q[5], a], s.pipe[5], s.fluid, s.conv)
    dQ[5] = pipe_cheng_valve_fit(t, [P[5], P[6], Q[5], a], s.pipe[5], s.fluid, s.conv)
    dQ[6] = pipe_cheng(t, [P[6], P_o, Q[6]], s.pipe[6], s.fluid, s.conv)
    # Pressure
    dP[1] = pipe_p(t, [Q[1], Q[2]], s.pipe[1], s.fluid, s.conv)
    dP[2] = pipe_p(t, [Q[2], Q[3]], s.pipe[2], s.fluid, s.conv)
    dP[3] = pipe_p(t, [Q[3], Q_p], s.pipe[3], s.fluid, s.conv)
    dP[4] = pipe_p(t, [Q_p, Q[4]], s.pipe[4], s.fluid, s.conv)
    dP[5] = pipe_p(t, [Q[4], Q[5]], s.pipe[5], s.fluid, s.conv)
    dP[6] = pipe_p(t, [Q[5], Q[6]], s.pipe[6], s.fluid, s.conv)
end


# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Twin pump
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1 - Turbulent - Cheng - Shock - Valve 

function hqvp1i1tcvfit!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    s = p.sys
    Q_p, ω = view(u, 1:2)
    dQ, Q = view(du, 3:4), view(u, 3:4)
    dP, P = view(du, 5:6), view(u, 5:6)
    # Impeller with laminar friction
    du[1] = impeller_fit1(t, [P[1], P[2], Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    du[2] = shaft_1st_fit(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    dQ[1] = pipe_cheng_twin_fit(t, [P_i, P[1], Q[1]], s.pipe[1], s.fluid, s.conv)
    # Pipeline downstream
    dQ[2] = pipe_cheng_valve_fit(t, [P[2], P_o, Q[2], a], s.pipe[2], s.fluid, s.conv)
    # Pressure
    dP[1] = pipe_p(t, [Q[1], Q_p], s.pipe[1], s.fluid, s.conv)
    dP[2] = pipe_p(t, [Q_p, Q[2]], s.pipe[2], s.fluid, s.conv)
end

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Twin pump
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1 - Laminar - Cheng - Shock - Valve 

function hqvp1i1lcvfit!(du, u, p, t)
    P_i, P_o = 0.0, 0.0
    s = p.sys
    Q_p, ω = view(u, 1:2)
    dQ, Q = view(du, 3:4), view(u, 3:4)
    dP, P = view(du, 5:6), view(u, 5:6)
    # Impeller with laminar friction
    du[1] = impeller_fit1(t, [P[1], P[2], Q_p, ω], s.pump, s.fluid, s.conv)
    # Shaft with first order friction only
    du[2] = shaft_1st_fit(t, [Q_p, ω], p.input, s.shaft, s.pump, s.fluid, s.conv)
    # Pipeline upstream
    dQ[1] = pipe_lam_twin(t, [P_i, P[1], Q[1]], s.pipe[1], s.fluid, s.conv)
    # Pipeline downstream
    dQ[2] = pipe_lam_valve(t, [P[2], P_o, Q[2], a], s.pipe[2], s.fluid, s.conv)
    # Pressure
    dP[1] = pipe_p(t, [Q[1], Q_p], s.pipe[1], s.fluid, s.conv)
    dP[2] = pipe_p(t, [Q_p, Q[2]], s.pipe[2], s.fluid, s.conv)
end