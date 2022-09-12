using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations
using ModelingToolkit

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

@named suc = Junction1(P)
@named pm = Junction0()

@named imp = Junction1(m)

@named val = Junction1([-1, d])
@named pj = Junction0()

cons = [connect(suc.power, pm.power), connect(lek.power, pm.power), connect(pm.power, imp.power), connect(imp.power, pj.power), connect(pj.power, lek.power), connect(pj.power, val.power)]

#cons = [connect(suc.power, pm.power), connect(pm.power, imp.power), connect(imp.power, pj.power), connect(pj.power, val.power)]
@named psys = ODESystem(cons, t)
mdl = compose(psys, lek, suc, pm, imp, pj, val)
#mdl = compose(psys, suc, pm, imp, pj, val)

generate_graph(mdl)

emdl = expand_connections(mdl)
equations(emdl)

sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)

latexify(equations(sys))

equations(sys)
