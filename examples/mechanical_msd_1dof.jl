# # [Harmonic Oscillator](@id msd-1dof)
#
# A mass, spring and damper obey `I u̇ + R u + q/C = F(t)` and `q̇ = u`.
# The source bond points into the 1-junction; the I, C and R bonds point out.
#
# ```@raw html
# <picture style="display: flex; justify-content: center; align-items: center;">
#   <source srcset="../../assets/tikz/bg_1dof-dark.svg" media="(prefers-color-scheme: dark)">
#   <img src="../../assets/tikz/bg_1dof.svg">
# </picture>
# ```

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t
using OrdinaryDiffEqTsit5
using Plots

function oscillator(force, displacement; name)
    @named m = Inertia(I = 1.0)
    @named s = Compliance(C = 1.0, q = displacement)
    @named d = Resistance(R = 1.0)
    @named f = Se(force)
    Junction1(f, [-1, m], [-1, d], [-1, s]; name)
end

# The unforced oscillator starts with unit spring displacement.
@named unforced = oscillator(0.0, 1.0)
unforced_sys = mtkcompile(unforced)
unknowns(unforced_sys)
equations(unforced_sys)
unforced_sol = solve(ODEProblem(unforced_sys, [], (0.0, 10.0)), Tsit5())
plot(unforced_sol, idxs = [unforced.s.q, unforced.m.power.f],
     xlabel = "Time", ylabel = "Displacement / velocity")

# A sinusoidal effort source drives the oscillator from rest.
@named forced = oscillator(sin(t), 0.0)
forced_sol = solve(ODEProblem(mtkcompile(forced), [], (0.0, 10.0)), Tsit5())
plot(forced_sol, idxs = [forced.s.q, forced.m.power.f],
     xlabel = "Time", ylabel = "Displacement / velocity")
