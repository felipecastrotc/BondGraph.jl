using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations
using ModelingToolkit

# Notes
# [15:09, 05/08/2022] Jorge Biazussi: Rampa de partida
# [15:09, 05/08/2022] Jorge Biazussi: 1hz/segundo
# [15:09, 05/08/2022] Jorge Biazussi: Até 40hz
# [15:10, 05/08/2022] Jorge Biazussi: Rampa de aceleração de até 3hz/seg

# ==========================================================
# Used alias
@variables ω(t), Q(t)

# ==========================================================
# New test with multiple outputs for pump leakage

# ----------------------------------------------------------
# Elements
@parameters p

@named m = Mass()
@named d = Damper()
@named s = Spring()
@named P = Se(p)

# ----------------------------------------------------------
# Junctions

@named lek = Junction1([-1, d])

@named ilt = Junction1(P)
@named jds = Junction0()

@named imp = Junction1(m)

@named olt = Junction1([-1, d])
@named jus = Junction0()

# ----------------------------------------------------------
# Connections

cons = [
    connect(ilt.power, jds.power),
    connect(lek.power, jds.power),
    connect(jds.power, imp.power),
    connect(imp.power, jus.power),
    connect(jus.power, lek.power),
    connect(jus.power, olt.power),
]

# ----------------------------------------------------------
# System

@named psys = ODESystem(cons, t)
mdl = compose(psys, lek, imp, jds, jus, ilt, olt)
#mdl = compose(psys, suc, pm, imp, pj, val)

generate_graph(mdl)

emdl = expand_connections(mdl)
equations(emdl)

sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)
equations(sys)

# ==========================================================
# Without pump leakage but with Jorge's losses

# ==========================================================
# Define losses expressions
@variables ω(t), Q(t)
@parameters ρ, d, μ, k1, k2, k3, k4, k5, k6, k7, k8, n
q = Q

# Friction losses -> Jorge eq 3.11 -> Look dev_pump_algebrics.jl
ΔPₐ₁ = k2 * q * μ / d^3
ΔPₐ₂ = k3 * q^(2 - n) * ρ^(1 - n) * d^(n - 4) * μ^n
# BG - Divided by Q for BG default R*Q
ΔPₐ₁ = ΔPₐ₁ / q
ΔPₐ₂ = k3 * q^(1 - n) * ρ^(1 - n) * d^(n - 4) * μ^n

# Vortex losses -> Jorge eq 3.15 -> Look dev_pump_algebrics.jl
# It is important to note that it possibly includes recirculation effects
# Look Paternost work
ΔPₜ₁ = d^2 * k4 * ρ * ω^2
ΔPₜ₂ = -2 * k4 * k5 * q * ρ * ω / d
ΔPₜ₃ = k4 * k5^2 * q^2 * ρ / d^4
# BG - Divided by Q for BG default R*Q
ΔPₜ₁ = ΔPₜ₁ / q
ΔPₜ₂ = ΔPₜ₂ / q
ΔPₜ₃ = ΔPₜ₃ / q

# Localized loss -> Jorge eq 3.16 -> Look dev_pump_algebrics.jl
ΔPₗ = k6 * q^2 * ρ / d^4
# BG - Divided by Q for BG default R*Q
ΔPₗ = ΔPₗ / q

# ----------------------------------------------------------
# Elements

# Losses
@named dΔPa1 = Damper(c = ΔPₐ₁)
@named dΔPa2 = Damper(c = ΔPₐ₂)
@named dΔPt1 = Damper(c = ΔPₜ₁)
@named dΔPt2 = Damper(c = ΔPₜ₂)
@named dΔPt3 = Damper(c = ΔPₜ₃)
@named dΔPl = Damper(c = ΔPₗ)

@named m = Mass()
@named s = Spring()
@named P = Se(p)

# ----------------------------------------------------------
# Junctions

@named lek = Junction1([-1, d])

@named ilt = Junction1(P)
@named jds = Junction0()

@named imp = Junction1(m)

@named olt = Junction1([-1, d])
@named jus = Junction0()

# ----------------------------------------------------------
# Connections

cons = [
    connect(ilt.power, jds.power),
    connect(lek.power, jds.power),
    connect(jds.power, imp.power),
    connect(imp.power, jus.power),
    connect(jus.power, lek.power),
    connect(jus.power, olt.power),
]

# ----------------------------------------------------------
# System

@named psys = ODESystem(cons, t)
mdl = compose(psys, lek, imp, jds, jus, ilt, olt)
#mdl = compose(psys, suc, pm, imp, pj, val)

generate_graph(mdl)

emdl = expand_connections(mdl)
equations(emdl)

sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)
equations(sys)

kg * m^2 / s^2


kg^2 / m
