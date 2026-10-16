# # [MSD system with N-DoF](@id msd-ndof)
#
# A chain of masses can be built by repeating the same bond pattern. Each mass has
# an inertia on a 1-junction; the first also has a compliance and resistance to
# ground. A compliance–resistance pair couples each pair of neighbouring masses.
# All parameters are one, and this example uses four degrees of freedom.
#
# A named free port represents one end of a bond. `[-1, :out]` points out of a mass
# junction; `:in` points into the next mass junction. At each plain `connect`, those
# ports and `sd.j0` form a 0-junction: effort is shared and signed flows sum to zero.
# The coupling bond points from that 0-junction into `sd` (`0 → sd`) because `:j0`
# is positive at `sd`. This is the textbook orientation. Relative to the v0.2
# drawing, that bond is reversed, which flips the spring displacement `q`. To keep
# the same initial stretch,
# we therefore start the second coupling compliance at `q = -1`.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t
using OrdinaryDiffEqTsit5
using Plots

n = 4
@named m = Inertia(I = 1.0)
@named s = Compliance(C = 1.0)
@named d = Resistance(R = 1.0)
@named m1 = Junction1([-1, m], [-1, d], [-1, s], [-1, :out])

masses = [Junction1(:in, [-1, m], (i < n ? [[-1, :out]] : [])...;
                    name = Symbol(:m, i)) for i in 2:n]
sds = [Junction1([-1, Compliance(C = 1.0, q = i == 2 ? -1.0 : 0.0, name = :s)],
                 [-1, d], :j0; name = Symbol(:sd, i)) for i in 2:n]
prev = [m1; masses[1:end-1]]
cons = [connect(prev[k].out, masses[k].in, sds[k].j0) for k in 1:n-1]
@named model = System(cons, t; systems = [m1; masses; sds])

sys = mtkcompile(model)
unknowns(sys)
sol = solve(ODEProblem(sys, [], (0.0, 10.0)), Tsit5())
plot(sol, xlabel = "Time", ylabel = "Amplitude")
