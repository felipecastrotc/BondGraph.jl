using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
# using ModelingToolkit
# using DifferentialEquations


# =============================================================================
# Mass dampser spring system
# -----------------------------------------------------------------------------
# Definitions

m = 1.0
X = 1.0
k = 10.0
c = 1.0
@named mass = Mass(m = m)
@named spring = Spring(k = k, x = X)
@named damper = Damper(c = c)

@named msd = Junction1(-mass, -spring, -damper, couple = false)

equations(alias_elimination(msd))


equations(msd)
sys = structural_simplify(msd)
@named syso = reducedobs(structural_simplify(msd))

equations(sys)
states(sys)
observed(sys)

equations(syso)
states(syso)
