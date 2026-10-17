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
@named mass = Mass(m = 2.0, u = 0.5)
@named spring = Spring(k = 2.0, x = 0.2)
@named damper = Damper(c = 3.0)

# A custom nonlinear compliance uses the same `power` port as the built-ins.
@oneport function CubicCompliance(; name, C = 1.0, q = 0.0)
    @parameters C = C
    @variables q(t) = q
    [power.e ~ q^3 / C, D(q) ~ power.f]
end
@named cubic = CubicCompliance(C = 1.0)

@named model = Junction1(effort_source, [-1, inertia], [-1, compliance],
                         [-1, resistance], [-1, cubic])
sol = solve(ODEProblem(mtkcompile(model), [], (0.0, 1.0)), Tsit5())
sol(1.0; idxs = model.compliance.q)

# ## Writing an element without the macros
#
# `@oneport` and `@twoport` only add the ports and wrap the laws. Any ModelingToolkit system with `PowerPort`s is an element, so a multiport, or an element with extra connectors, can be written by hand. This nonlinear resistance, a turbulent orifice with `e = R f |f|`, is written in its conductance form `f = sign(e) √(|e|/R)`, so the flow is explicit in the effort. It is the same system `@oneport` would build:
@component function QuadraticResistance(; name, R = 1.0)
    @named power = PowerPort()
    @parameters R = R
    System([power.f ~ sign(power.e) * sqrt(abs(power.e) / R)], t; name, systems = [power])
end
@named orifice = QuadraticResistance(R = 2.0)
@named se = Se(8.0)
@named model = Junction1(se, [-1, orifice])
sol = solve(ODEProblem(mtkcompile(model), [], (0.0, 1.0)), Tsit5())
sol(1.0; idxs = model.orifice.power.f)   # f = √(E/R) = 2
