# # [Reusable Subsystems](@id composition)
#
# A subsystem written once and used several times is a Julia function that returns a model with a power port. This page builds a series RC load, uses it twice, reads results from each copy, and changes one copy's parameter without recompiling.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t
using OrdinaryDiffEqTsit5
using Plots

# ## A one-junction subsystem
#
# A junction is already a subsystem: it owns its elements, and a free port is its connection to the outside. The series RC load is a 1-junction whose free port `power` is a bond entering it (`:power => :in`), so the load is connected like an element, through `load.power`:
RCLoad(; name, R = 1.0, C = 1.0) =
    Junction1(:power => :in, Resistance(; R, name = :r), Compliance(; C, name = :c); name)

# Two copies, each driven by its own source. Every copy gets its own namespace, `A.r`, `B.r`, so their elements, states and parameters stay apart:
@named A = RCLoad(R = 2.0, C = 0.5)
@named B = RCLoad(R = 3.0, C = 0.25)
@named VA = Se(1.0)
@named VB = Se(2.0)
@named loads = System([connect(VA.pout, A.power), connect(VB.pout, B.power)], t;
                      systems = [VA, VB, A, B])
bgplot(loads)

# ## Results from each copy
#
# Results are read through the model's names, as for any element: `loads.A.c.q` is the charge of A's capacitor.
sys = mtkcompile(loads)
prob = ODEProblem(sys, [], (0.0, 5.0))
sol = solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-10)
plot(sol, idxs = [loads.A.c.q, loads.B.c.q], labels = ["q, A" "q, B"], xlabel = "Time (s)")

# By hand, a source $V$ charging $R$ and $C$ in series from empty gives $q(t) = C V (1 - e^{-t/RC})$. The largest difference on a common time grid, for A and for B:
q_rc(τ; R, C, V) = C * V * (1 - exp(-τ / (R * C)))
tgrid = 0:0.01:5
(maximum(abs, sol(tgrid; idxs = loads.A.c.q).u - q_rc.(tgrid; R = 2.0, C = 0.5, V = 1.0)),
 maximum(abs, sol(tgrid; idxs = loads.B.c.q).u - q_rc.(tgrid; R = 3.0, C = 0.25, V = 2.0)))

# ## Changing one copy's parameter
#
# `remake` gives the problem new parameter values without compiling again. Name the parameter through the **compiled** system, `sys.A.r.R`. The model before compiling, `loads`, still carries its own name in its paths (`loads₊A₊r₊R`), which solution indexing accepts but `remake` rejects with a warning that it "is not a parameter". This is ModelingToolkit's behaviour: `complete(loads).A.r.R` works too.
prob_slow = remake(prob; p = [sys.A.r.R => 4.0])

# `prob.ps` reads a parameter's value. Only A's resistance changed:
(prob.ps[sys.A.r.R], prob_slow.ps[sys.A.r.R], prob_slow.ps[sys.B.r.R])

# A charges more slowly now, and B is unchanged:
sol_slow = solve(prob_slow, Tsit5(), reltol = 1e-10, abstol = 1e-10)
plot(sol, idxs = [loads.A.c.q], label = "A, R = 2", xlabel = "Time (s)", ylabel = "Charge")
plot!(sol_slow, idxs = [loads.A.c.q, loads.B.c.q], labels = ["A, R = 4" "B"], linestyle = :dash)

# Against the closed form with A's new resistance, and B against its first run:
(maximum(abs, sol_slow(tgrid; idxs = loads.A.c.q).u - q_rc.(tgrid; R = 4.0, C = 0.5, V = 1.0)),
 maximum(abs, sol_slow(tgrid; idxs = loads.B.c.q).u - sol(tgrid; idxs = loads.B.c.q).u))

# ## A subsystem with several junctions
#
# A subsystem with more than one junction is a plain ModelingToolkit component that owns a `PowerPort` and connects it inside, as here an R on a 1-junction followed by a C on a 0-junction. The port's direction comes from the inner port it is connected to: `j1.power` is a bond entering the 1-junction, so `ladder.power` is one too.
@component function RCLadder(; name, R = 1.0, C = 1.0)
    @named power = PowerPort()
    @named j1 = Junction1(:power => :in, Resistance(; R, name = :r), :out)
    @named j0 = Junction0(:in => :in, Compliance(; C, name = :c))
    System([connect(power, j1.power), connect(j1.out, j0.in)], t; name, systems = [power, j1, j0])
end
@named ladder = RCLadder(R = 2.0, C = 0.5)
@named V = Se(1.0)
@named driven = System([connect(V.pout, ladder.power)], t; systems = [V, ladder])
sol = solve(ODEProblem(mtkcompile(driven), [], (0.0, 5.0)), Tsit5(), reltol = 1e-10, abstol = 1e-10)
maximum(abs, sol(tgrid; idxs = driven.ladder.j0.c.q).u - q_rc.(tgrid; R = 2.0, C = 0.5, V = 1.0))

# The drawing shows such a component as one node, named after it: `bondgraph` doesn't draw the junctions and elements inside a component it doesn't know. A one-junction subsystem, like `RCLoad` above, is drawn in full, because it is a junction.
bgplot(driven)
