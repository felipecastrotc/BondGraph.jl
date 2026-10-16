# # [Duffing Equation](@id msd-duffing)
#
# The Duffing oscillator adds a cubic compliance to the linear mass, compliance
# and resistance: `I u̇ + R u + q/C + q³/C₃ = F(t)`. Both compliances share one
# 1-junction. The source bond points into it; all element bonds point out.
#
# ```@raw html
# <picture style="display: flex; justify-content: center; align-items: center;">
#   <source srcset="../../assets/tikz/bg_duffing-dark.svg" media="(prefers-color-scheme: dark)">
#   <img src="../../assets/tikz/bg_duffing.svg">
# </picture>
# ```

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using Plots

# Custom one-ports use the bond variables `power.e` and `power.f` directly.
@oneport function CubicCompliance(; name, C = 0.2, q = 0.0)
    @parameters C = C
    @variables q(t) = q
    [power.e ~ q^3 / C, D(q) ~ power.f]
end

function duffing(force, displacement; name)
    @named m = Inertia(I = 1.0)
    @named s = Compliance(C = 1.0, q = displacement)
    @named s3 = CubicCompliance(C = 0.2, q = displacement)
    @named d = Resistance(R = 0.02)
    @named f = Se(force)
    Junction1(f, [-1, m], [-1, d], [-1, s], [-1, s3]; name)
end

@named unforced = duffing(0.0, 1.0)
unforced_sys = mtkcompile(unforced)
unknowns(unforced_sys)
equations(unforced_sys)
unforced_sol = solve(ODEProblem(unforced_sys, [], (0.0, 10.0)), Tsit5())
plot(unforced_sol, idxs = [unforced.s.q, unforced.m.power.f],
     xlabel = "Time", ylabel = "Displacement / velocity")

@named forced = duffing(8.0cos(0.5t), 0.0)
forced_sol = solve(ODEProblem(mtkcompile(forced), [], (0.0, 10.0)), Tsit5())
plot(forced_sol, idxs = [forced.s.q, forced.m.power.f],
     xlabel = "Time", ylabel = "Displacement / velocity")
