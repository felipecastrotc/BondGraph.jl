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

# -------------------------------------------------------------------
# Pipe models
# -------------------------------------------------------------------
# Pipeline with laminar friction only
function pipe_lam(t, u, pipe, fluid, conv)
    # Unpack states
    P_in, P_out, Q, Q_out = u[:]
    # Unpack properties
    μ = fluid.μ
    l, d, I, C = pipe.l, pipe.d, pipe.I, pipe.C
    # Unpack useful constants
    H = conv.H
    # Equation
    return [
        (P_in - P_out - darcyq_lam(l, d, μ) * Q / H) * H / I,
        (Q - Q_out) / (C * H)
    ]
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