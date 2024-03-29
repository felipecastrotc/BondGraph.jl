using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D

# ============================================================================
# Mass-spring-damper problems

# Definitions
@named m = Mass(m = 1)
@named s = Spring(k = 1, x = 0.1)
@named s3 = Spring3(k = 1, x = 0.2)
@named d = Damper(c = 1)
@named F = Se(5)

# 1-DOF
# @named dof1 = Junction1(m, s, d, F)
@named dof1 = Junction1([-1, m], [-1, s], [-1, d], F)
generate_graph(dof1)
edof1 = expand_connections(dof1)

equations(edof1)

@named sdof1 = reducedobs(structural_simplify(edof1))
equations(sdof1)
prob = ODEProblem(sdof1, [0.0, 0.0], (0.0, 10))
sol = solve(prob)
plot(sol)

# 2-DOF
@named sd = Junction1(s, d, F)
@named mj = Junction1(m)

@named b1 = Junction0(mj, sd)
@named b0 = Junction1(m, s, d)

@named dof2 = ODESystem([connect(b0.power, b1.power)], t)
cdof2 = compose(dof2, b1, b0)
generate_graph(cdof2)
edof2 = expand_connections(cdof2)

@named sdof2 = reducedobs(structural_simplify(edof2))
equations(sdof2)

states(sdof2)
prob = ODEProblem(sdof2, [0.0, 0.0, 0.0, 0.0], (0.0, 10))
sol = solve(prob)
plot(sol)

# n-DOF
@named sd = Junction1(s, d)
@named mj = Junction1(m)

@named b1 = Junction0(mj, sd)
@named b2 = Junction0(mj, sd)

@named b0 = Junction1(m, s, d)

@named dofn = ODESystem([connect(b2.power, b1.mj.power), connect(b1.power, b0.power)], t)
cdofn = compose(dofn, b1, b2, b0)
generate_graph(cdofn)

edofn = expand_connections(cdofn)
@named sdofn = reducedobs(structural_simplify(edofn))
equations(sdofn)

# ============================================================================
# MGY problems

@named L = Mass(m = 0.5)
@named R = Damper(c = 1.0)
@named Uₐ = Se(-12.0)

@named J = Mass(m = 0.01)
@named b = Damper(c = 0.1)
@named Tₗ = Se(-1.0)

g = 0.01

@named jm = Junction1(Tₗ, b, J)
@named je = Junction1(Uₐ, R, L)
# @named jm = Junction1(b, J)
@named gy = mGY(je, jm, g = g)
equations(gy)

@named mdl = compose(gy, je, jm)
generate_graph(mdl)
emdl = expand_connections(mdl)
@named sys = reducedobs(structural_simplify(emdl))
equations(sys)

# Example coherent with literature
# https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling
