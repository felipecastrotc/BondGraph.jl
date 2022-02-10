using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D

using DifferentialEquations


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

