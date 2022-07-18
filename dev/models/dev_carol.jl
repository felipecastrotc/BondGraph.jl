using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations, JLD2


@named se = Se(10);
@named im = Mass(m=1.0);
@named rm = Damper(c=2.0);
@named c = Spring(k=3.0);
@named i2 = Mass(m=4.0);
@named ro = Damper(c=5.0);

@named j11 = Junction1(se, -rm, -im);
@named j0 = Junction0(j11, -c);
@named sys = Junction0(j0, -i2, -ro, couple=false);

@named mdl = reducedobs(structural_simplify(sys));

equations(mdl)