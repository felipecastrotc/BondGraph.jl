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

# ----------------------------------------------------------
# Parameters and properties
@parameters Q0, P0

a = 1500                # m/s water
ρ = 998                 # kg/mˆ3 water
μ = 1e-3                # Pa*s water
a = 1290                # m/s petroleum
ρ = 850                 # kg/mˆ3 petroleum
μ = 800e-3              # Pa*s petroleum
d = 75 / 1000             # m Pipe diameter
Δs = 7.5                # m Pipe length
ϵ = 0.1 * 1e-3            # m Steel rugosity
A = π * (d / 2)^2       # m^2 Pipe area

q = 1e-4                # Flowrate
q = 100e-4                # Flowrate
Re = (ρ * (q / A) * d) / μ
ReEq(ρ, q, d, μ)

L = ρ * Δs / A
C = A * Δs / (ρ * a^2)
R = 32 * μ / (A * d^2) * Δs
# Δs * cheng(Re, ϵ, d) * (ρ / (2 * d))/(A^2)
# R = 5e7

# Check these values, including Reynolds with
# http://www.druckverlust.de/Online-Rechner/dp.php
cheng(Re, ϵ, d)
colebrook(Re, ϵ, d)
Δs * cheng(Re, ϵ, d) * (ρ / (2 * d)) * ((q / A)^2)
darcyq(q, Δs, d, ϵ, ρ, μ) * (q^2)

# Bond graph elements
@named m = Mass(m = L)       # General inertial element
@named f = Damper(c = R)     # General damping element
@named s = Spring(k = -1 / C)   # General compliance element
@named Qi = Sf(Q0)           # Inlet pressure
@named Po = Se(P0)          # Inlet pressure

# ----------------------------------------------------------------
# Pipeline 1 element
@named s11 = Junction1([-1, f], [-1, m])
@named s10 = Junction0([-1, s], Qi)

cons = [connect(s10.power, s11.power)]
@named psys = ODESystem(cons, t)

mdl = compose(psys, s11, s10)
emdl = expand_connections(mdl)
sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)

equations(sys)
states(sys)
parameters(sys)

# The resistance is seriously problematic use the manual drop loss
vals = [s10.Qi.Q0 => 1]
v0 = [q, 0]
prob = ODEProblem(sys, v0, (0.0, 2), vals)

# sol = solve(prob, reltol=1e-8, abstol=1e-8)
# sol = solve(prob)
# plot(sol)

# ----------------------------------------------------------------
# Pipeline 2-element
@named s21 = Junction1([-1, f], [-1, m])
@named s20 = Junction0([-1, s])

push!(cons, connect(s11.power, s20.power), connect(s20.power, s21.power))
@named psys = ODESystem(cons, t)

mdl = compose(psys, s11, s10, s21, s20)
emdl = expand_connections(mdl)
sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)

generate_graph(mdl)
equations(sys)
states(sys)
parameters(sys)

# The resistance is seriously problematic use the manual drop loss
vals = [s10.Qi.Q0 => 1]
v0 = [q, 0]
prob = ODEProblem(sys, v0, (0.0, 2), vals)

# sol = solve(prob, reltol=1e-8, abstol=1e-8)
# sol = solve(prob)
# plot(sol)

# ----------------------------------------------------------------
# Pipeline 3-element
@named s31 = Junction1([-1, f], [-1, m])
@named s30 = Junction0([-1, s])

push!(cons, connect(s21.power, s30.power), connect(s30.power, s31.power))
@named psys = ODESystem(cons, t)

mdl = compose(psys, s11, s10, s21, s20, s30, s31)
emdl = expand_connections(mdl)
sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)

generate_graph(mdl)
equations(sys)
states(sys)
parameters(sys)

# The resistance is seriously problematic use the manual drop loss
vals = [s10.Qi.Q0 => 1]
v0 = [q, 0]
prob = ODEProblem(sys, v0, (0.0, 2), vals)

# sol = solve(prob, reltol=1e-8, abstol=1e-8)
# sol = solve(prob)
# plot(sol)

# --------------------------------------------------------------
# FORGET BOND-GRAPH
# --------------------------------------------------------------
# From this time I give up from modifying the sytem and modelig manually
# Bond graph is still useful on how to couple things

# For lower reynolds number it might be interesting to increse manually
# increase the resistance by 10 and reduce the Compliance (C) by
# multiplying it by 1000

Kc = 50
Kf = 2

# Manual 1-element
function pipe!(du, u, p, t)
    Q2, P1 = u
    du[1] = ((P1 - 0) - Kf * darcyq(Q2, Δs, d, ϵ, ρ, μ) * (Q2^2)) / L
    du[2] = (q - Q2) / (Kc * C)
end

q
u0 = [q * 1.5; 0.0]
tspan = (0.0, 2)
prob = ODEProblem(pipe!, u0, tspan)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

plot(sol)
plot(sol.t, sol[1, :])
plot(sol.t, sol[2, :])

# Manual 2-element
function pipe!(du, u, p, t)
    Q2, Q4, P1, P3 = u
    du[1] = (P1 - P3 - Kf * darcyq(Q2, Δs, d, ϵ, ρ, μ) * (Q2^2)) / L
    du[2] = (P3 - Kf * darcyq(Q4, Δs, d, ϵ, ρ, μ) * (Q4^2)) / L
    du[3] = (q - Q2) / (Kc * C)
    du[4] = (Q2 - Q4) / (Kc * C)
end

u0 = [q * 1.5; q * 1.5; 0.0; 0.0]
tspan = (0.0, 2)
prob = ODEProblem(pipe!, u0, tspan)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

plot(sol)
plot(sol.t, sol[1, :])
plot!(sol.t, sol[2, :])
plot(sol.t, sol[3, :])
plot!(sol.t, sol[4, :])

# Manual 3-element
function pipe!(du, u, p, t)
    Q2, Q4, Q6, P1, P3, P5 = u
    du[1] = (P1 - P3 - Kf * darcyq(Q2, Δs, d, ϵ, ρ, μ) * (Q2^2)) / L
    du[2] = (P3 - P5 - Kf * darcyq(Q4, Δs, d, ϵ, ρ, μ) * (Q4^2)) / L
    du[3] = (P5 - Kf * darcyq(Q6, Δs, d, ϵ, ρ, μ) * (Q6^2)) / L
    du[4] = (q - Q2) / (Kc * C)
    du[5] = (Q2 - Q4) / (Kc * C)
    du[6] = (Q4 - Q6) / (Kc * C)
end

u0 = [q * 1.5; q * 1.5; q * 1.5; 0.0; 0.0; 0.0]
tspan = (0.0, 3)
prob = ODEProblem(pipe!, u0, tspan)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

plot(sol)
plot(sol.t, sol[1, :])
plot!(sol.t, sol[2, :])
plot!(sol.t, sol[3, :])

plot(sol.t, sol[4, :])
plot!(sol.t, sol[5, :])
plot!(sol.t, sol[6, :])
