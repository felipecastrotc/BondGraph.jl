# # [Mass-Spring-Damper system 2-DoF](@id msd-2dof)
#
# Two masses exchange force through a coupling spring and resistance. The first mass
# also has a spring and resistance to ground; the second receives an external force.
# We use unit inertias, compliances, and resistances, so each compliance represents a
# spring with stiffness $k = 1/C = 1$.
#
# Each 1-junction represents a shared velocity. A named free port exposes one bond:
# `[-1, :out]` points out of the first mass junction, while `:in` points into the
# second mass and coupling spring–resistance junctions. Their plain `connect` is the
# intervening 0-junction: it shares effort and balances the signed flows. The signs
# belong to the far ends of the bonds, not to the order of `connect` arguments.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t
using OrdinaryDiffEqTsit5
using Plots

# The unforced case starts with unit displacement in the coupling spring. The forced
# case starts at rest and applies $F(t)=\sin(t)$ to the second mass.
for (label, force, q) in (("Unforced", 0.0, 1.0), ("Forced", sin(t), 0.0))
    @named m = Inertia(I = 1.0)
    @named s = Compliance(C = 1.0)
    @named d = Resistance(R = 1.0)
    @named f = Se(force)
    @named m1 = Junction1([-1, m], [-1, d], [-1, s], [-1, :out])
    @named m2 = Junction1([-1, m], f, :in)
    @named s2 = Compliance(C = 1.0, q = q)
    @named sd = Junction1([-1, s2], [-1, d], :in)
    @named model = System([connect(m1.out, m2.in, sd.in)], t; systems = [m1, m2, sd])

    sys = mtkcompile(model)
    unknowns(sys)
    sol = solve(ODEProblem(sys, [], (0.0, 10.0)), Tsit5())
    plot(sol, title = label, xlabel = "Time", ylabel = "Amplitude")
end
