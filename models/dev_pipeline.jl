using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations
using ModelingToolkit

# ==========================================================
# Used alias
@variables ω(t), Q(t)

# ----------------------------------------------------------
# Elements
@parameters p, po, τ, Qₚ, ρ

Ap = π * (3 * 0.0254 / 2)^2
m1 = 30 * 998 / (Ap)

@named m = Mass(m = rpo["Iₚᵢ"])       # General inertial element
@named m = Mass(m = m1)       # General inertial element
@named d1 = Damper(c = rpo["fₚᵢ"] * 20)     # General damping element
# @named d1 = Damper(c=0.2)     # General damping element
@named d2 = Damper(c = rpo["fₚᵢ"] * 20)     # General damping element
@named d3 = Damper(c = rpo["fₚᵢ"] * 20)     # General damping element
@named s1 = Spring(k = -1 / rpo["kᵤ"])     # General compliance element
# @named s1 = Spring(k=-0.1)     # General compliance element
@named s2 = Spring(k = -1 / rpo["kᵤ"])     # General compliance element
@named s3 = Spring(k = -1 / rpo["kᵤ"])     # General compliance element
@named P = Se(p)        # Inlet pressure
@named Po = Se(p)        # Inlet pressure
@named P = Se(p)        # Inlet pressure
@named Po = Se(p)        # Inlet pressure

# Pipeline
@named s11 = Junction1([-1, d1], [-1, m], P)
@named s10 = Junction0([-1, s1])

@named s21 = Junction1([-1, d2], [-1, m])
@named s20 = Junction0([-1, s2])

@named s31 = Junction1([-1, d3], [-1, m])
@named s30 = Junction0([-1, s3], [-1, Po])

cons = [
    connect(s11.power, s10.power),
    connect(s10.power, s21.power),
    connect(s21.power, s20.power),
    connect(s20.power, s31.power),
    connect(s31.power, s30.power),
]
# cons = [connect(s11.power, s10.power), connect(s10.power, s21.power), connect(s21.power, s20.power)]
# cons = [connect(s11.power, s10.power)]

@named psys = ODESystem(cons, t)

# mdl = compose(psys, s11, s10, s21, s20, s31)
mdl = compose(psys, s11, s10, s21, s20, s31, s30)
# mdl = compose(psys, s11, s10, s21, s20)
# mdl = compose(psys, s11, s10)

# mdl = s11
generate_graph(mdl)

emdl = expand_connections(mdl)
equations(emdl)
amdl = ModelingToolkit.alias_elimination(emdl)
equations(amdl)
sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)

equations(sys)
states(sys)
parameters(sys)

vals = [s11.P.p => 1e3]
vals = [s11.P.p => 4e3]
vals = [s11.P.p => 4e3, s30.Po.p => 0.0]
# vals = [s11.P.p => 4e3, s30.Po.p => 4e3]
# vals = [P.p => 4e3]
v0 = [0, 0, 0, 0, 0]
# v0 = [0, 0, 0]
# v0 = [0, 0]
# v0 = [0]
prob = ODEProblem(sys, v0, (0.0, 30), vals)
equations(sys)

sol = solve(prob, reltol = 1e-13, abstol = 1e-13)
# sol = solve(prob)
plot(sol)
sol ./ Ap
equations(sys)

plot(sol.t, sol[s11.m.power.f])
plot(sol.t, sol[s10.s1.q])
