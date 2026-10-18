# # [One-Port Elements](@id elements-oneport)
#
# A one-port exchanges effort and flow along its bond. Sources conventionally
# point away from themselves; passive elements conventionally point into them.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5

@named effort_source = Se(1.0)
@named flow_source = Sf(1.0)
@named inertia = Inertia(I = 2.0, p = 1.0)
@named compliance = Compliance(C = 0.5, q = 0.2)
@named resistance = Resistance(R = 3.0)

# Mechanical shorthands use mass and initial velocity, stiffness and initial
# displacement, and damping coefficient respectively.
@named mass = Mass(m = 2.0, v = 0.5)
@named spring = Spring(k = 2.0, x = 0.2)
@named damper = Damper(c = 3.0)

# A custom nonlinear compliance gets the same ports (`power`, `pout`) and bond variables (`e`, `f`) as the built-ins. These four names are taken inside `@oneport`: call a frequency `ν`, not `f`.
@oneport function CubicCompliance(; name, C = 1.0, q = 0.0)
    @parameters C = C
    @variables q(t) = q
    [e ~ q^3 / C, D(q) ~ f]
end
@named cubic = CubicCompliance(C = 1.0)

@named model = Junction1(effort_source => :pout, inertia, compliance,
                         resistance, cubic)
sol = solve(ODEProblem(mtkcompile(model), [], (0.0, 1.0)), Tsit5())
sol(1.0; idxs = model.compliance.q)

# ## Writing an element without the macros
#
# `@oneport` and `@twoport` only add the ports and wrap the laws. Any ModelingToolkit system with `PowerPort`s is an element, so a multiport, or an element with extra connectors, can be written by hand. This nonlinear resistance, a turbulent orifice with `e = R f |f|`, is written in its conductance form `f = sign(e) √(|e|/R)`, so the flow is explicit in the effort. It is the same system `@oneport` would build: the incoming port `power`, the outgoing port `pout`, and the bond's `e` and `f`, whichever port carries the bond (ModelingToolkit zeroes the flow of the port left unconnected):
@component function QuadraticResistance(; name, R = 1.0)
    @named power = PowerPort()
    @named pout = PowerPort()
    @variables e(t) f(t)
    @parameters R = R
    System([e ~ power.e, pout.e ~ power.e, f ~ power.f - pout.f,
            f ~ sign(e) * sqrt(abs(e) / R)], t; name, systems = [power, pout])
end
@named orifice = QuadraticResistance(R = 2.0)
@named se = Se(8.0)
@named model = Junction1(se => :pout, orifice)
sol = solve(ODEProblem(mtkcompile(model), [], (0.0, 1.0)), Tsit5())
sol(1.0; idxs = model.orifice.f)   # f = √(E/R) = 2
