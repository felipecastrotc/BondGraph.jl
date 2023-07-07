using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations

@named se = Se(10);
@named im = Mass(m=1.0);
@named rm = Damper(c=2.0);
@named c = Spring(k=3.0);
@named i2 = Mass(m=4.0);
@named ro = Damper(c=5.0);

@named j1 = Junction1(se, [-1, rm], [-1, im]);
@named j01 = Junction0([-1, c]);
@named j02 = Junction0([-1, i2], [-1, ro]);

eqs = [connect(j1.power, j01.power), connect(j01.power, j02.power)]
@named sys = ODESystem(eqs, t)
sys = compose(sys, j1, j01, j02)

@named mdl = reducedobs(structural_simplify(sys));

equations(mdl)
generate_graph(sys, method=:sfdp)