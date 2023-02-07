using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations, JLD2


# Pag 383 Karnop

# ----------------------------------------------------------
# Parameters and properties
@parameters Q0, P0, P2

# Bond graph elements
@named m = Mass()       # General inertial element
@named f = Damper()     # General damping element
@named s = Spring()   # General compliance element
@named ω = Sf(Q0)           # Inlet pressure
@named Pi = Se(P0)          # Inlet pressure
@named Po = Se(P2)          # Inlet pressure


# Upstream
@named up0 = Junction0(Pi)
# Downstream
@named dw0 = Junction0()
# Leakage loss
@named lk = Junction1([-1, f])
con = [connect(up0.power, lk.power), connect(lk.power, dw0.power)]

# Pump parameters
@parameters k1, k2

@named c1 = mTF(r=k1)
@named c2 = mTF(r=k1)

@named p1 = Junction1(ω)

push!(con, connect(up0.power, c1.pin), connect(c1.pout, p1.power), connect(p1.power, c2.pin), connect(c2.pout, dw0.power))

# Pipe
@named s11 = Junction1([-1, f], [-1, m])
@named s10 = Junction0([-1, s], Po)

push!(con, connect(dw0.power, s11.power), connect(s11.power, s10.power))

# Define the system
@named psys = ODESystem(con, t)
# mdl = compose(psys, up0, dw0, lk, c1, c2, p1)
mdl = compose(psys, up0, dw0, lk, c1, c2, p1, s11, s10)
emdl = expand_connections(mdl)

equations(ModelingToolkit.alias_elimination(emdl))
sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)

generate_graph(mdl)

equations(sys)

# 1-element Vector{Equation}:
#  Differential(t)(s11₊m₊power₊f(t)) ~ (up0₊Pi₊P0 + lk₊f₊R*(-s11₊m₊power₊f(t) - k1*p1₊ω₊Q0) - s10₊Po₊P2 - s11₊f₊R*s11₊m₊power₊f(t)) / s11₊m₊I

