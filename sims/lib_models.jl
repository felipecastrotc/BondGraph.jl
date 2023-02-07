# -------------------------------------------------------------------
# Impeller models
# -------------------------------------------------------------------
# Impeller with laminar friction only
function impeller_lam(t, u, pump, fluid, conv)
    # Unpack states
    P_in, P_out, Q, ω = u[:]
    # Unpack properties
    ρ, μ = fluid.ρ, fluid.μ
    γa, γb, d, l, I = pump.γa, pump.γb, pump.d, pump.l, pump.I
    # Unpack useful constants
    H = conv.H
    # Equation
    return (P_in - P_out + (ρ / H) * (γa * Q - γb * ω) * ω - darcyq_lam(l, d, μ) * Q / H) * H / I
end

# Impeller with laminar friction only
function impeller_lshock(t, u, pump, fluid, conv)
    # Unpack states
    P_in, P_out, Q, ω = u[:]
    # Unpack properties
    ρ, μ = fluid.ρ, fluid.μ
    γa, γb, d, l, I = pump.γa, pump.γb, pump.d, pump.l, pump.I
    ks, q_star = pump.ks, pump.q_star
    # Unpack useful constants
    H = conv.H
    # Equation
    return (P_in - P_out + (ρ / H) * (γa * Q - γb * ω) * ω - darcyq_lam(l, d, μ) * Q / H - ks * ρ * (d / 2)^2 * ω^2 * (1 - Q / q_star)^2 / H) * H / I
end

# Impeller fitted - expanded with shock-loss expanded
function impeller_fit1(t, u, pump, fluid, conv)
    # Unpack states
    P_in, P_out, Q, ω = u[:]
    # Unpack properties
    ρ, μ = fluid.ρ, fluid.μ
    I = pump.I
    k_1, k_2, k_3, k_4 = pump.K
    c_f = pump.c_f
    # Unpack useful constants
    H = conv.H
    # Equation
    return (P_in - P_out + (k_1 * ρ * ω * Q + k_2 * ρ * ω^2 + k_3 * μ * Q + k_4 * ρ * (Q^2) + c_f) / H) * H / I
end

# -------------------------------------------------------------------
# Pipe models
# -------------------------------------------------------------------
# Pipeline with laminar friction only
function pipe_lam(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, I, k = pipe.l, pipe.d, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Unpack useful constants
    H = conv.H
    # Equation
    return (P_in - P_out - darcyq_lam(l, d, μ) * Q / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end

# Pipeline with haaland friction
function pipe_haal(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, ϵ, I, k = pipe.l, pipe.d, pipe.ϵ, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Unpack useful constants
    H = conv.H
    # Equation
    return (P_in - P_out - darcyq_haal(Q, l, d, ϵ, μ, ρ) * Q^2 / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end

# Pipeline with cheng friction
function pipe_cheng(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, ϵ, I, k = pipe.l, pipe.d, pipe.ϵ, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Unpack useful constants
    H = conv.H
    # Equation
    return (P_in - P_out - darcyq_cheng(Q, l, d, ϵ, μ, ρ) * Q^2 / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end

# Pipeline with cheng friction with valve
function pipe_cheng_valve(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q, a = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, ϵ, I, k = pipe.l, pipe.d, pipe.ϵ, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Valve properties
    C_v_func, C_v_ps = pipe.others.C_v, pipe.others.C_v_ps
    N_1, F_p = pipe.others.N_1, pipe.others.F_p
    # Unpack useful constants
    H = conv.H
    # Valve equations
    C_v = C_v_func(a, C_v_ps)
    ΔPv = ΔPValve(Q * 3600, ρ, C_v, N_1, F_p) * 1e3     # Conversions m^3/h and kPa
    # Pipe equations
    return (P_in - P_out - darcyq_cheng(Q, l, d, ϵ, μ, ρ) * Q^2 / H - ΔPv / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end

function pipe_lam_valve(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q, a = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, ϵ, I, k = pipe.l, pipe.d, pipe.ϵ, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Valve properties
    C_v_func, C_v_ps = pipe.others.C_v, pipe.others.C_v_ps
    N_1, F_p = pipe.others.N_1, pipe.others.F_p
    # Unpack useful constants
    H = conv.H
    # Valve equations
    C_v = C_v_func(a, C_v_ps)
    ΔPv = ΔPValve(Q * 3600, ρ, C_v, N_1, F_p) * 1e3     # Conversions m^3/h and kPa
    # Pipe equations
    return (P_in - P_out - darcyq_lam(l, d, μ) * Q / H - ΔPv / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end

# Pipeline with cheng friction with valve
function pipe_cheng_valve_fit(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q, a = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, ϵ, I, k = pipe.l, pipe.d, pipe.ϵ, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Valve properties
    C_v_func, C_v_ps = pipe.others.C_v, pipe.others.C_v_ps
    N_1, F_p = pipe.others.N_1, pipe.others.F_p
    # Unpack useful constants
    H = conv.H
    # Valve equations
    C_v = C_v_func(a, C_v_ps)
    ΔPv = ΔPValve(Q * 3600, ρ, C_v, N_1, F_p) * 1e3     # Conversions m^3/h and kPa
    # Pipe equations
    return (P_in - P_out - darcyq_cheng(Q, l, d, ϵ, μ, ρ) * Q^2 / H - ΔPv / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end

# Pipeline with cheng friction with cavity pump
function pipe_cheng_cavity(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, ϵ, I, k = pipe.l, pipe.d, pipe.ϵ, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Cavity pump properties
    k1, k2, ω = pipe.others.k_1, pipe.others.k_2, pipe.others.ω
    # Unpack useful constants
    H = conv.H
    # Cavity pump equation
    P_0 = pcavity(Q, μ, ω, k1, k2, 0)
    # Pipe equations
    return (P_0 / H - P_out - darcyq_cheng(Q, l, d, ϵ, μ, ρ) * Q^2 / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end

# Pipeline with cheng friction with twin-screw pump
function pipe_cheng_twin(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, ϵ, I, k = pipe.l, pipe.d, pipe.ϵ, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Twin-screw pump properties
    k1, k2, ω, h = pipe.others.k_1, pipe.others.k_2, pipe.others.ω, pipe.others.h
    # Unpack useful constants
    H = conv.H
    # Twin-screw pump equation
    ωi = smooth_startup_10(t, ω)
    P_0 = ptwin(Q, μ, ωi, k1, k2, ρ * 9.81 * h)
    # Pipe equations
    return (P_0 / H - P_out - darcyq_cheng(Q, l, d, ϵ, μ, ρ) * Q^2 / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end

function pipe_lam_twin(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, ϵ, I, k = pipe.l, pipe.d, pipe.ϵ, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Twin-screw pump properties
    k1, k2, ω, h = pipe.others.k_1, pipe.others.k_2, pipe.others.ω, pipe.others.h
    # Unpack useful constants
    H = conv.H
    # Twin-screw pump equation
    ωi = smooth_startup_10(t, ω)
    P_0 = ptwin(Q, μ, ωi, k1, k2, ρ * 9.81 * h)
    # Pipe equations
    return (P_0 / H - P_out - darcyq_lam(l, d, μ) * Q / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end

# Pipeline with cheng friction with twin-screw pump fit
function pipe_cheng_twin_fit(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q = u[:]
    # Unpack properties
    μ, ρ = fluid.μ, fluid.ρ
    l, d, ϵ, I, k = pipe.l, pipe.d, pipe.ϵ, pipe.I, pipe.k
    A = π * ((d / 2) ^ 2)
    # Twin-screw pump properties
    k1, k2, ω, h = pipe.others.k_1, pipe.others.k_2, pipe.others.ω, pipe.others.h
    # Unpack useful constants
    H = conv.H
    # Twin-screw pump equation
    P_0 = ptwin(Q, μ, ω, k1, k2, ρ * 9.81 * h)
    # Pipe equations
    return (P_0/H - P_out - darcyq_cheng(Q, l, d, ϵ, μ, ρ) * Q^2 / H - (k * (ρ/2)*(Q/A)^2)/H) * H / I
    # return (P_0 / H - P_out - darcyq_lam(l, d, μ) * Q / H - (k * (ρ/2)*(Q/A)^2) / H) * H / I
end


# Pipeline pressure
function pipe_p(t, u, pipe, fluid, conv)
    # Unpack states
    Q_in, Q_out = u[:]
    # Unpack properties
    C = pipe.C
    # Unpack useful constants
    H = conv.H
    # Equation
    return (Q_in - Q_out) / (C * H)
end

# -------------------------------------------------------------------
# Shaft model
# -------------------------------------------------------------------
# Shaft with first order friction only
function shaft_1st(t, u, input, shaft, pump, fluid, conv)
    # Unpack states
    Q, ω = u[:]
    # Unpack properties
    ρ = fluid.ρ
    γa, γb, = pump.γa, pump.γb
    c, I = shaft.c, shaft.I
    # Unpack input function
    τ, τp = input.func, input.p
    # Unpack useful constants
    H = conv.H
    # Equation
    return (τ(t, τp...) - c * ω - ρ * (γa * Q - γb * ω) * Q) / I
end

# Shaft with first order friction only
function shaft_1st_fit(t, u, input, shaft, pump, fluid, conv)
    # Unpack states
    Q, ω = u[:]
    # Unpack properties
    ρ, μ = fluid.ρ, fluid.μ
    I = shaft.I
    k_1, k_2, k_3, k_4, k_5 = shaft.K
    # Unpack input function
    τ, τp = input.func, input.p
    # Unpack useful constants
    H = conv.H
    # Equation
    return (τ(t, τp...) - k_1*ρ*(Q^2) - k_2*ρ*ω*Q - k_3*μ*ω - k_4*ω - k_5*(ω^2)) / I
end