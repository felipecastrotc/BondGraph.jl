# # [Coupling with the MTK Standard Library](@id adapters)
#
# ModelingToolkitStandardLibrary (MTKSL) models use their own connectors: pins, flanges, hydraulic ports. A bond graph can't `connect` to them directly, because their variables are not an effort and a flow. An *adapter* is a one-port bond-graph element whose other side is an MTKSL connector, so the MTKSL network behind it acts as the element at the end of the bond:
#
# | Adapter | MTKSL connector(s) | Effort | Flow |
# |---|---|---|---|
# | `ElectricalAdapter` | `p`, `n` (`Pin`) | voltage `p.v - n.v` | current into the network at `p` |
# | `TranslationalAdapter` | `flange` (`MechanicalPort`) | force on the network | `flange.v` |
# | `TranslationalPositionAdapter` | `flange` (`TranslationalPosition.Flange`) | force on the network | `ds/dt` |
# | `RotationalAdapter` | `flange` (`Rotational.Flange`) | torque on the network | `dφ/dt` |
# | `HydraulicAdapter` | `port` (`HydraulicPort`) | pressure | volume flow into the network |
#
# Like any element, an adapter's law is written in the variables of its bond, and its bond conventionally points into the adapter: `[-1, adapter]` at a junction. When the MTKSL side supplies power instead of absorbing it, the power along the bond is negative.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t
import ModelingToolkitStandardLibrary.Electrical as El
import ModelingToolkitStandardLibrary.Mechanical.Rotational as Rot
import ModelingToolkitStandardLibrary.Mechanical.Translational as Tr
import ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition as TP
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as Hy
import ModelingToolkitStandardLibrary.Blocks as Blocks
using OrdinaryDiffEqRosenbrock
using Plots

# ## A bond-graph motor driving an MTKSL load
#
# The [DC motor](@ref dc-motor) with its rotor inertia and friction moved to MTKSL: the mechanical 1-junction bonds into a `RotationalAdapter`, whose flange drives an MTKSL `Inertia` and `Damper`. The adapter integrates the shaft speed into the flange angle, so it takes the initial angle `phi`.
@named L = Inertia(I = 0.5)
@named R = Resistance(R = 1.0)
@named V = Se(12.0)
@named τ = Se(1.0)
@named shaft = RotationalAdapter(phi = 0.0)
@named je = Junction1(V, [-1, R], [-1, L], [-1, :gy])
@named jm = Junction1(τ, :gy, [-1, shaft])
@named gy = GY(r = 0.01)
@named inertia = Rot.Inertia(J = 0.01, w = 0.0)
@named damper = Rot.Damper(d = 0.1)
@named fixed = Rot.Fixed()
@named motor = System([connect(je.gy, gy.p1), connect(gy.p2, jm.gy),
                       connect(jm.shaft.flange, inertia.flange_a), connect(inertia.flange_b, damper.flange_b),
                       connect(damper.flange_a, fixed.flange)], t;
                      systems = [je, gy, jm, inertia, damper, fixed])
sol = solve(ODEProblem(mtkcompile(motor), [], (0.0, 10.0)), Rodas5P(), reltol = 1e-10, abstol = 1e-10);

# As a reference, the same motor written entirely in MTKSL: a voltage source, resistor and inductor drive an `EMF`, which couples them like the gyrator (torque `k i`, back EMF `k ω`), whose flange turns the same inertia and damper, with a constant applied torque. The `EMF` integrates the shaft speed into an angle, so here the rotor takes its initial angle `phi`.
@named vk = Blocks.Constant(k = 12.0)
@named src = El.Voltage()
@named res = El.Resistor(R = 1.0)
@named ind = El.Inductor(L = 0.5, i = 0.0)
@named emf = El.EMF(k = 0.01)
@named rotor = Rot.Inertia(J = 0.01, phi = 0.0, w = 0.0)
@named gnd = El.Ground()
@named τk = Blocks.Constant(k = 1.0)
@named torque = Rot.Torque()
@named ref = System([connect(vk.output, src.V), connect(src.p, res.p), connect(res.n, ind.p),
                     connect(ind.n, emf.p), connect(emf.n, src.n, gnd.g),
                     connect(τk.output, torque.tau), connect(emf.flange, rotor.flange_a, torque.flange),
                     connect(rotor.flange_b, damper.flange_b), connect(damper.flange_a, emf.support, fixed.flange)], t;
                    systems = [vk, src, res, ind, emf, gnd, τk, torque, rotor, damper, fixed])
sol_ref = solve(ODEProblem(mtkcompile(ref), [], (0.0, 10.0)), Rodas5P(), reltol = 1e-10, abstol = 1e-10)
plot(sol, idxs = [inertia.w, je.L.power.f], label = ["ω, bond graph + MTKSL" "i, bond graph + MTKSL"])
plot!(sol_ref, idxs = [rotor.w, ind.i], linestyle = :dash, label = ["ω, MTKSL" "i, MTKSL"], xlabel = "Time")

# The largest difference on a common time grid, in angular velocity and in current:
tgrid = 0:0.01:10
(maximum(abs, sol(tgrid; idxs = inertia.w).u - sol_ref(tgrid; idxs = rotor.w).u),
 maximum(abs, sol(tgrid; idxs = je.L.power.f).u - sol_ref(tgrid; idxs = ind.i).u))

# ## An MTKSL circuit behind a bond graph source
#
# Here a bond-graph source drives an MTKSL resistor and capacitor in series through an `ElectricalAdapter`; the network's lower terminal is grounded.
@named se = Se(3.0)
@named ea = ElectricalAdapter()
@named j = Junction1(se, [-1, ea])
@named res = El.Resistor(R = 2.0)
@named cap = El.Capacitor(C = 0.5, v = 0.0)
@named gnd = El.Ground()
@named circuit = System([connect(j.ea.p, res.p), connect(res.n, cap.p), connect(cap.n, j.ea.n, gnd.g)], t;
                        systems = [j, res, cap, gnd])
sol = solve(ODEProblem(mtkcompile(circuit), [], (0.0, 5.0)), Rodas5P(), reltol = 1e-10, abstol = 1e-10)

# The reference is the same circuit with an MTKSL voltage source in place of the bond-graph source and adapter:
@named vk = Blocks.Constant(k = 3.0)
@named src = El.Voltage()
@named ref = System([connect(vk.output, src.V), connect(src.p, res.p), connect(res.n, cap.p),
                     connect(cap.n, src.n, gnd.g)], t; systems = [vk, src, res, cap, gnd])
sol_ref = solve(ODEProblem(mtkcompile(ref), [], (0.0, 5.0)), Rodas5P(), reltol = 1e-10, abstol = 1e-10)
plot(sol, idxs = [cap.v], label = "Bond graph + MTKSL", xlabel = "Time", ylabel = "Capacitor voltage")
plot!(sol_ref, idxs = [cap.v], linestyle = :dash, label = "MTKSL")

# The largest difference on a common time grid:
tgrid = 0:0.01:5
maximum(abs, sol(tgrid; idxs = cap.v).u - sol_ref(tgrid; idxs = cap.v).u)

# ## Translational: a force on an MTKSL mass
#
# A bond-graph force pushes a 2 kg MTKSL mass, first velocity-based (`Translational`), then position-based (`TranslationalPosition`). By Newton's 2nd law, a 1 N force accelerates the mass from rest as `v = t/2`, and both models are checked against it. The velocity-based mass needs its own initial position and velocity:
@named F = Se(1.0)
@named ta = TranslationalAdapter()
@named j = Junction1(F, [-1, ta])
@named mass = Tr.Mass(m = 2.0)
@named push = System([connect(j.ta.flange, mass.flange)], t; systems = [j, mass])
sys = mtkcompile(push)
sol = solve(ODEProblem(sys, [sys.mass.v => 0.0, sys.mass.s => 0.0], (0.0, 2.0)), Rodas5P(),
            reltol = 1e-10, abstol = 1e-10)
tgrid = 0:0.01:2
maximum(abs, sol(tgrid; idxs = mass.v).u - tgrid / 2)
#-
@named F = Se(1.0)
@named tpa = TranslationalPositionAdapter()
@named j = Junction1(F, [-1, tpa])
@named mass = TP.Mass(m = 2.0, s = 0.0, v = 0.0)
@named push = System([connect(j.tpa.flange, mass.flange)], t; systems = [j, mass])
sol = solve(ODEProblem(mtkcompile(push), [], (0.0, 2.0)), Rodas5P(), reltol = 1e-10, abstol = 1e-10)
maximum(abs, sol(tgrid; idxs = mass.v).u - tgrid / 2)

# ## Hydraulic: an MTKSL pressure source filling a bond-graph tank
#
# Here the MTKSL side supplies power: a pressure source drives a flow through a bond-graph resistance into a compliance (a tank). With the adapter's bond drawn into it, the bond's flow is negative while the source delivers, so the tank's volume `q` goes negative: the volume that left the source.
@named k = Blocks.Constant(k = 2.0e5)
@named src = Hy.Pressure()
@named fluid = Hy.HydraulicFluid()
@named ha = HydraulicAdapter()
@named R = Resistance(R = 1.0e9)
@named C = Compliance(C = 2.0e-9)
@named j = Junction1([-1, ha], [-1, R], [-1, C])
@named tank = System([connect(k.output, src.p), connect(fluid, src.port), connect(src.port, j.ha.port)], t;
                     systems = [k, src, fluid, j])
sol = solve(ODEProblem(mtkcompile(tank), [], (0.0, 10.0)), Rodas5P(); abstol = 1e-18, reltol = 1e-10)

# By hand: all three bonds point out of the 1-junction, so its efforts sum to $0 = -p_0 - R\,\dot{q} - q/C$, and from an empty tank
#
# $q(t) = -C\,p_0\,\left(1 - e^{-t/RC}\right),$
#
# which approaches $-C\,p_0 = -4 \times 10^{-4}$ m³ with the time constant $RC = 2$ s.
q_tank(τ) = -2.0e-9 * 2.0e5 * (1 - exp(-τ / (1.0e9 * 2.0e-9)))
plot(sol, idxs = [j.C.q], label = "Bond graph + MTKSL", xlabel = "Time", ylabel = "Volume (m³)")
plot!(q_tank, 0, 10, label = "Closed form", linestyle = :dash)

# The largest difference on a common time grid, relative to the final volume $C\,p_0$:
tgrid = 0:0.01:10
maximum(abs, sol(tgrid; idxs = j.C.q).u - q_tank.(tgrid)) / 4.0e-4

# Notes per domain:
#
# - `TranslationalPositionAdapter` and `RotationalAdapter` integrate the flow into a position or angle, so give their initial value (`s`, `phi`) unless the MTKSL side already sets it.
# - `HydraulicAdapter`: MTKSL's through-variable is a mass flow, so the bond's flow is the volume flow at the fluid's reference density. Connect the `HydraulicFluid` in its own statement, `connect(fluid, src.port)`, not in the same `connect` as the adapter. Hydraulic volumes are small (~1e-6 m³), so tighten `abstol` accordingly.
# - A coupled model may keep algebraic equations. If a solver reports that it can't use mass matrices, switch to a DAE-capable one such as `Rodas5P`, as here.
