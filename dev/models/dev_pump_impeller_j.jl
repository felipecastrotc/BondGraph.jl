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
# Define losses expressions
@variables ω(t), Q(t)
@parameters ρ, d, μ, k1, k2, k3, k4, k5, k6, k7, k8, n

# Friction losses -> Jorge eq 3.11 -> Look dev_pump_algebrics.jl
k2*q*μ/d^3 + k3*q^2*ρ*(d*μ/(q*ρ))^n/d^4
ΔPₐ₁ = k2*q*μ/d^3
ΔPₐ₂ = k3*q^(2-n)*ρ^(1-n)*d^(n-4)*μ^n

# Vortex losses -> Jorge eq 3.15 -> Look dev_pump_algebrics.jl
# It is important to note that it possibly includes recirculation effects
# Look Paternost work
ΔPₜ₁ = d^2*k4*ρ*ω^2
ΔPₜ₂ = - 2*k4*k5*q*ρ*ω/d
ΔPₜ₃ = k4*k5^2*q^2*ρ/d^4

# Localized loss -> Jorge eq 3.16
# ΔPₗ = k6*ρ*(ω^2)*(d^2)*(Q/(ω*d^3))^2
ΔPₗ = k6*ρ*(ω^2)*(d^2)*(Q^2/(ω^2*d^6))    # simplified form from above


# ==========================================================
# Elements
@parameters p

@named m = Mass()
@named d = Damper()
@named s = Spring()
@named P = Se(p)

# ==========================================================
# New test with multiple outputs for pump leakage

@named lek = Junction1([-1, d])

@named ilt = Junction1(P)
@named jds = Junction0()

@named imp = Junction1(m)

@named olt = Junction1([-1, d])
@named jus = Junction0()

cons = [connect(ilt.power, jds.power), connect(lek.power, jds.power), connect(jds.power, imp.power), connect(imp.power, jus.power), connect(jus.power, lek.power), connect(jus.power, olt.power)]

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
# Adding Jorge thesis model losses

exp(0.2*log(0))

0^(0.5)
