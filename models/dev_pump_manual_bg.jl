using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations
using ModelingToolkit

function ReEq(ρ, q, d, μ)
    A = π * (d / 2)^2
    return (ρ * (q / A) * d) / μ
end

function darcyq(q, Δs, d, ϵ, ρ, μ)
    Re = ReEq(ρ, q, d, μ)
    if Re > 0.0
        return Δs * cheng(Re, ϵ, d) * (ρ / (2 * d)) * ((1 / A)^2)
    else
        return 0.0
    end
end


function colebrook(Re, ϵ, d)
    # Return the friction factor as 0 for reynolds number equals to 0
    if Re < 2300
        return min(64 / Re, 64)
    end
    # Niazkar approximation
    # https://link.springer.com/article/10.1007/s12205-019-2217-1
    A = -2 * log10((ϵ / d) / 3.7 + 4.5547 / (Re^0.8784))
    B = -2 * log10((ϵ / d) / 3.7 + 2.51 * A / (Re))
    C = -2 * log10((ϵ / d) / 3.7 + 2.51 * B / (Re))
    # Estimate f
    rhs = A - ((B - A)^2) / (C - 2 * B + A)
    f = 1 / (rhs^2)

    # # Set the Niazkar approximation as an initial guess to iterate
    # f₀ = f
    # for i in 1:20
    #     # Colebrook equation
    #     f = (1 / (-2 * log10(ϵ / (3.7 * d) + 2.51 / (Re * sqrt(f₀)))))^2
    #     # Absolute error with tolerance of 1e-7
    #     if abs(f - f₀) < 1e-7
    #         continue
    #     end
    #     f₀ = f
    # end
    return f
end

# Friction factor
function cheng(Re, ϵ, d)
    a = 1 / (1 + (Re / 2720)^9)
    b = 1 / (1 + (Re / (160 * (d / ϵ)))^2)
    invf =
        ((Re / 64)^a) *
        ((1.8 * log10(Re / 6.8))^(2 * (1 - a) * b)) *
        ((2.0 * log10(3.7 * d / ϵ))^(2 * (1 - a) * (1 - b)))
    return 1 / invf
end


# ==========================================================
# Used alias
@variables ω(t), Q(t)

# ------------------------------------------------------------------------------
# Elements
@parameters Q0, P0, τ

@named m = Mass()       # General inertial element
@named r = Damper()     # General damping element
@named s = Spring()     # General compliance element
@named Qi = Sf(Q0)           # Inlet pressure
@named Po = Se(P0)          # Inlet pressure
@named Pi = Se(P0)          # Inlet pressure
@named T = Se(τ)        # Shaft input torque

# ------------------------------------------------------------------------------
# Impeller Junctions

@named lek = Junction1([-1, r])     # Leakage setting

@named jus = Junction0()     # Coupling for leakage

@named imp = Junction1([-1, m], [-1, r])     # Impeller

@named olt = Junction1([-1, r])     # Shock losses
@named jds = Junction0()            # Coupling for leakage

# ------------------------------------------------------------------------------
# Axis Junctions
@named ax = Junction1(T, [-1, m], [-1, r])      # Axis parameters

# ------------------------------------------------------------------------------
# Pipe Junctions -  source of flow

# # Inlet
# @named us11 = Junction1([-1, r], [-1, m])
# @named us10 = Junction0([-1, s], Qi)

# @named us21 = Junction1([-1, r], [-1, m])
# @named us20 = Junction0([-1, s])

# plist = [us11, us10, us21, us20]

# @named ds11 = Junction1([-1, r], [-1, m])
# @named ds10 = Junction0([-1, s])

# @named ds21 = Junction1([-1, r], [-1, m], Po)
# @named ds20 = Junction0([-1, s])

# push!(plist, ds11, ds10, ds21, ds20);

# ------------------------------------------------------------------------------
# Pipe Junctions -  source of effort

# Inlet
@named us11 = Junction1([-1, r], [-1, m], Pi)
@named us10 = Junction0([-1, s])

@named us21 = Junction1([-1, r], [-1, m])
@named us20 = Junction0([-1, s])

@named us31 = Junction1([-1, r], [-1, m])
@named us30 = Junction0([-1, s])

@named us41 = Junction1([-1, r], [-1, m])
@named us40 = Junction0([-1, s])

# plist = [us11, us10, us21, us20]
plist = [us11, us10, us21, us20, us31, us30, us41, us40]

@named ds11 = Junction1([-1, r], [-1, m])
@named ds10 = Junction0([-1, s])

@named ds21 = Junction1([-1, r], [-1, m])
# @named ds20 = Junction0([-1, s], Po)
@named ds20 = Junction0([-1, s])

@named ds31 = Junction1([-1, r], [-1, m])
@named ds30 = Junction0([-1, s])

@named ds41 = Junction1([-1, r], [-1, m])
@named ds40 = Junction0([-1, s], Po)

push!(plist, ds11, ds10, ds21, ds20);
push!(plist, ds31, ds30, ds41, ds40);


# ------------------------------------------------------------------------------
# MGY
@parameters γa, γb, ρ

gyp = ρ * (γa * ω - γb * Q)
@named agy = mGY(g = gyp)

# ------------------------------------------------------------------------------
# Build connections

# Pipeline assembly - source of flow
# cons = [connect(us10.power, us11.power), connect(us11.power, us20.power), connect(us20.power, us21.power)]
# push!(cons, connect(ds10.power, ds11.power), connect(ds11.power, ds20.power), connect(ds20.power, ds21.power))
# Pipeline assembly - source of effort
cons = [
    connect(us11.power, us10.power),
    connect(us10.power, us21.power),
    connect(us21.power, us20.power),
]
push!(
    cons,
    connect(ds11.power, ds10.power),
    connect(ds10.power, ds21.power),
    connect(ds21.power, ds20.power),
)
# 4 elements
push!(
    cons,
    connect(us20.power, us31.power),
    connect(us31.power, us30.power),
    connect(us30.power, us41.power),
    connect(us41.power, us40.power),
)
push!(
    cons,
    connect(ds20.power, ds31.power),
    connect(ds31.power, ds30.power),
    connect(ds30.power, ds41.power),
    connect(ds41.power, ds40.power),
)

# Impeller assembly - source of flow
# push!(cons, connect(us21.power, jus.power), connect(jus.power, imp.power), connect(imp.power, jds.power), connect(jds.power, olt.power), connect(olt.power, ds10.power))
# Impeller assembly - source of effort
# push!(cons, connect(us20.power, jus.power), connect(jus.power, imp.power), connect(imp.power, jds.power), connect(jds.power, olt.power), connect(olt.power, ds10.power))
# 4 elements pipe
push!(
    cons,
    connect(us40.power, jus.power),
    connect(jus.power, imp.power),
    connect(imp.power, jds.power),
    connect(jds.power, olt.power),
    connect(olt.power, ds10.power),
)

# Axis - Impeller coupling
push!(cons, connect(ax.power, agy.pin), connect(agy.pout, imp.power))
push!(cons, agy.Q ~ agy.pin.f, agy.ω ~ agy.pout.f)

# ------------------------------------------------------------------------------
# System

@named psys = ODESystem(cons, t)
mdl = compose(psys, imp, jds, jus, olt, ax, agy, plist...)
# mdl = compose(psys, imp, jds, jus, olt, plist...)

generate_graph(mdl)

emdl = expand_connections(mdl)
equations(emdl)
sys2 = ModelingToolkit.structural_simplify(emdl)
@named sys2 = reducedobs(sys2)
equations(sys2)
equations(sys)
# --------------------------------------------------------------
# FORGET BOND-GRAPH
# --------------------------------------------------------------
# From this time I give up from modifying the sytem and modelig manually
# Bond graph is still useful on how to couple things

# For lower reynolds number it might be interesting to increse manually
# increase the resistance by 10 and reduce the Compliance (C) by
# multiplying it by 1000

# Properties
# Fluid
# a = 1500                # m/s water
# ρ = 998                 # kg/mˆ3 water
# μ = 1e-3                # Pa*s water
a = 1290                # m/s petroleum
ρ = 850                 # kg/mˆ3 petroleum
μ = 800e-3              # Pa*s petroleum
# Pipeline
d = 75 / 1000             # m Pipe diameter
Δs = 7.5                # m Pipe length
ϵ = 0.1 * 1e-3            # m Steel rugosity
A = π * (d / 2)^2       # m^2 Pipe area
# Pipeline dynamics -> using Tanaka pipe
# https://sci-hub.se/10.1109/IECON.2000.972506
L = ρ * Δs / A
C = A * Δs / (ρ * a^2)
R = 32 * μ / (A * d^2) * Δs

# State conditions
q = 1e-4
# Reynolds check
Re = (ρ * (q / A) * d) / μ
ReEq(ρ, q, d, μ)

# Friction check
# Check these values, including Reynolds with
# http://www.druckverlust.de/Online-Rechner/dp.php
cheng(Re, ϵ, d)
colebrook(Re, ϵ, d)
Δs * cheng(Re, ϵ, d) * (ρ / (2 * d)) * ((q / A)^2)
darcyq(q, Δs, d, ϵ, ρ, μ) * (q^2)

# -------------------------------------------------------------------
# Model manual -  3 upstream 2 downstream - Source of flow
function pump_pipe!(du, u, p, t)
    Qi, ω, Q2, Q4, Q6, Q8, Q10, P1, P3, P5, P7, P9 = u[:]
    Δs, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L = p[:]
    # du = zeros(size(u))
    # Impeller
    du[1] =
        (P5 - P7 + ρ * (γa * Qi - γb * ω) * ω - darcyq(Qi, Δs, d, ϵ, ρ, μ) * (Qi^2)) / Ii
    # Shaft
    du[2] = (τ - ca * ω^2 - ρ * (γa * Qi - γb * ω) * ω) / Ia
    # Pipeline
    # Upstream
    # Q2
    du[3] = (P1 - P3 - Kf * darcyq(Q2, Δs, d, ϵ, ρ, μ) * (Q2^2)) / L
    # Q4
    du[4] = (P3 - P5 - Kf * darcyq(Q4, Δs, d, ϵ, ρ, μ) * (Q4^2)) / L
    # Q6
    du[5] = (P5 - Kf * darcyq(Q6, Δs, d, ϵ, ρ, μ) * (Q6^2)) / L
    # Downstream
    # Q8
    du[6] = (P7 - P9 - Kf * darcyq(Q8, Δs, d, ϵ, ρ, μ) * (Q8^2)) / L
    # Q10
    du[7] = (P9 - Po - Kf * darcyq(Q10, Δs, d, ϵ, ρ, μ) * (Q10^2)) / L
    # Pressures
    du[8] = (q - Q2) / (Kc * C)           # P1 Upstream inlet
    du[9] = (Q2 - Q4) / (Kc * C)          # P3 Upstream coupling pump
    du[10] = (Q4 - Qi) / (Kc * C)         # P5 Upstream coupling pump
    du[11] = (Qi - Q8) / (Kc * C)         # P9 Downstream coupling pump
    du[12] = (Q8 - Q10) / (Kc * C)         # P11 PDownstream outlet
end

# P100
γa = 0.002
γb = -1.7317
γa = 0.02
γb = -0.7317
Ii = 0.8504
Ia = 5e-4
ca = 0.001
Kf = 2
Kc = 50

# Inputs
τ = 110  # Torque
q = 0    # flow-rate
q = 1e-4    # flow-rate
Po = 0.0

q0 = q * 1.5
u0 = [0, 0, q0, q0, q0, q0, q0, 0, 0, 0, 0, 0]
p = [Δs, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L]

tspan = (0.0, 2.5)
prob = ODEProblem(pump_pipe!, u0, tspan, p)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

u0 = sol[:, end]

plot(sol.t, sol[1, :])
plot(sol.t, sol[2, :])
plot(sol.t, sol[3, :])
plot!(sol.t, sol[4, :])
plot!(sol.t, sol[5, :])

# -------------------------------------------------------------------
# Model manual -  2 upstream 2 downstream - Source of effort
function pump_pipe!(du, u, p, t)
    Qi, ω, Q1, Q3, Q5, Q7, P2, P4, P6, P8 = u[:]
    Δs, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L = p[:]
    # du = zeros(size(u))
    # Impeller
    du[1] =
        (P4 - P6 + ρ * (γa * Qi - γb * ω) * ω - darcyq(Qi, Δs, d, ϵ, ρ, μ) * (Qi^2)) / Ii
    # du[1] = (P2 - P6 + 8*ρ * (γa * Qi - γb * ω) * ω - darcyq(Qi, Δs, d, ϵ, ρ, μ) * (Qi^2)) / Ii
    # Shaft
    du[2] = (τ(t) - ca * ω - ρ * (γa * Qi - γb * ω) * Qi) / Ia
    # Pipeline
    # Upstream
    # Q1
    du[3] = (Pi - P2 - Kf * darcyq(Q1, Δs, d, ϵ, ρ, μ) * (Q1^2)) / L
    # Q3
    du[4] = (P2 - P4 - Kf * darcyq(Q3, Δs, d, ϵ, ρ, μ) * (Q3^2)) / L
    # Downstream
    # Q5
    du[5] = (P6 - Kf * darcyq(Q5, Δs, d, ϵ, ρ, μ) * (Q5^2)) / L
    # Q7
    du[6] = (P8 - Po - Kf * darcyq(Q7, Δs, d, ϵ, ρ, μ) * (Q7^2)) / L
    # Pressures
    du[7] = (Q1 - Q3) / (Kc * C)           # P2 Upstream inlet
    # du[7] = (Q1 - Qi) / (Kc * C)           # P2 Upstream inlet
    du[8] = (Q3 - Qi) / (Kc * C)           # P4 Upstream coupling pump
    du[9] = (Qi - Q5) / (Kc * C)           # P6 Downstream coupling pump
    du[10] = (Q5 - Q7) / (Kc * C)          # P8 PDownstream outlet
end

equations(sys2)
# P100
γa = 0.002
γb = -1.7317
# γa = 0.2
# γb = -10.7317
Ii = 0.8504
Ia = 5e-4
ca = 0.1
Kf = 2
Kc = 50

# Inputs
τ = (t) -> t > 0.5 ? 80 : 0.01    # Torque
# q = 0                      # flow-rate
# q = 1e-4                  # flow-rate
Po = 0
Pi = 0

# q = 1e-5
q = 0
q0 = q * 1.5
ω0 = 0.0
u0 = [0, ω0, q0, q0, q0, q0, 0, 0, 0, 0]
p = [Δs, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L]

tspan = (0.0, 2.5)
prob = ODEProblem(pump_pipe!, u0, tspan, p)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

60 * 12 / (2 * π)

# u0 = sol[:, end]
sol
# Q, ω
plot(sol.t, sol[1, :])
plot(sol.t, sol[2, :])
# Q1, Q3
plot(sol.t, sol[3, :])
plot!(sol.t, sol[4, :])
# Q5, Q7
plot(sol.t, sol[5, :])
plot!(sol.t, sol[6, :])
# P2, P4
plot(sol.t, sol[7, :])
plot!(sol.t, sol[8, :])
# P6, P8
plot(sol.t, sol[9, :])
plot!(sol.t, sol[10, :])

Qs = [1, 3, 4, 5, 6]
tl_q = [
    "Flow rate impeller",
    "Flow rate pipe segment 1",
    "Flow rate pipe segment 2",
    "Flow rate pipe segment 3",
    "Flow rate pipe segment 4",
]

Ps = [7, 8, 9, 10]
tl_p = [
    "Pressure drop pipe 1",
    "Pressure drop pipe 2",
    "Pressure drop pipe 3",
    "Pressure drop pipe 3",
    "Pressure drop pipe 4",
]

cfg = Dict(
    :titlefontsize => 8,
    :size => (800, 140 * length(Qs)),
    :legend => false,
    :xlabel => "Time (s)",
    :xguidefontsize => 8,
    :ylabel => "Q (m^3/s)",
    :yguidefontsize => 8,
)
plot(sol.t, sol[Qs, :]', layout = (length(Qs), 1), titles = reshape(tl_q, 1, 5); cfg...)

cfg = Dict(
    :titlefontsize => 8,
    :size => (800, 140 * length(Qs)),
    :legend => false,
    :xlabel => "Time (s)",
    :xguidefontsize => 8,
    :ylabel => "Pressure (Pa)",
    :yguidefontsize => 8,
)
plot(sol.t, sol[Ps, :]', layout = (length(Ps), 1), titles = reshape(tl_p, 1, 5); cfg...)

cfg = Dict(
    :titlefontsize => 8,
    :legend => false,
    :xlabel => "Time (s)",
    :xguidefontsize => 8,
    :ylabel => "Rotation (rad/s)",
    :yguidefontsize => 8,
)
plot(sol.t, sol[2, :], titles = "Shaft rotation velocity"; cfg...)

# -------------------------------------------------------------------
# Model manual -  4 upstream 4 downstream - Source of effort
function pump_pipe4!(du, u, p, t)
    Qi, ω, Q1, Q3, Q5, Q7, Q9, Q11, Q13, Q15, P2, P4, P6, P8, P10, P12, P14, P16 = u[:]
    Δs, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L = p[:]
    # du = zeros(size(u))
    # Impeller
    du[1] =
        (P8 - P10 + ρ * (γa * Qi - γb * ω) * ω - darcyq(Qi, Δs, d, ϵ, ρ, μ) * (Qi^2)) / Ii
    # Shaft
    du[2] = (τ(t) - ca * ω - ρ * (γa * Qi - γb * ω) * Qi) / Ia
    # Pipeline
    # Upstream
    # Q1
    du[3] = (Pi - P2 - Kf * darcyq(Q1, Δs, d, ϵ, ρ, μ) * (Q1^2)) / L
    # Q3
    du[4] = (P2 - P4 - Kf * darcyq(Q3, Δs, d, ϵ, ρ, μ) * (Q3^2)) / L
    # Q5
    du[5] = (P4 - P6 - Kf * darcyq(Q5, Δs, d, ϵ, ρ, μ) * (Q5^2)) / L
    # Q7
    du[6] = (P6 - P8 - Kf * darcyq(Q7, Δs, d, ϵ, ρ, μ) * (Q7^2)) / L
    # Downstream
    # Q9
    du[7] = (P10 - Kf * darcyq(Q9, Δs, d, ϵ, ρ, μ) * (Q9^2)) / L
    # Q11
    du[8] = (P12 - P14 - Kf * darcyq(Q11, Δs, d, ϵ, ρ, μ) * (Q11^2)) / L
    # Q13
    du[9] = (P14 - P16 - Kf * darcyq(Q13, Δs, d, ϵ, ρ, μ) * (Q13^2)) / L
    # Q15
    du[10] = (P16 - Po - Kf * darcyq(Q15, Δs, d, ϵ, ρ, μ) * (Q15^2)) / L
    # Pressures
    du[11] = (Q1 - Q3) / (Kc * C)           # P2  Upstream inlet
    du[12] = (Q3 - Q5) / (Kc * C)           # P4  Upstream inlet
    du[13] = (Q5 - Q7) / (Kc * C)           # P6  Upstream inlet
    du[14] = (Q7 - Qi) / (Kc * C)           # P8  Upstream coupling pump
    du[15] = (Qi - Q9) / (Kc * C)           # P10 Downstream coupling pump
    du[16] = (Q9 - Q11) / (Kc * C)          # P12 PDownstream outlet
    du[17] = (Q11 - Q13) / (Kc * C)         # P14 PDownstream outlet
    du[18] = (Q13 - Q15) / (Kc * C)         # P16 PDownstream outlet
end

equations(sys2)
# P100
γa = 0.002
γb = -1.7317
# γa = 0.2
# γb = -10.7317
Ii = 0.8504
Ia = 5e-4
ca = 0.1
Kf = 2
Kc = 5

# Inputs
τ = (t) -> t > 0.5 ? 2 : 0.0001    # Torque
# q = 0                      # flow-rate
# q = 1e-4                  # flow-rate
Po = 0
Pi = 0

# q = 1e-5
q = 0
q0 = q * 1.5
ω0 = 0.0
u0 = [0, ω0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
p = [Δs / 2, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L]

tspan = (0.0, 2.0)
prob = ODEProblem(pump_pipe4!, u0, tspan, p)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

plot(sol)

# Q, ω
plot(sol.t, sol[1, :])
plot(sol.t, sol[2, :])
# Q1, Q3, Q5, Q7
plot(sol.t, sol[3, :])
plot!(sol.t, sol[4, :])
plot!(sol.t, sol[5, :])
plot!(sol.t, sol[6, :])
# Q9, Q11, Q13
plot(sol.t, sol[7, :])
plot!(sol.t, sol[8, :])
plot!(sol.t, sol[9, :])
plot!(sol.t, sol[10, :])
