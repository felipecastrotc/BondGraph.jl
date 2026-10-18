# # [Equations and LaTeX](@id equations-latex)
#
# A bond-graph model is an MTK system, so its equations are available as for any MTK model. `bg_equations` returns them in terms of the states only, with the names you choose, ready for a report or a paper through `latexify`.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5

# A mass–spring–damper driven by a force: the source's bond points into the 1-junction, the element bonds point out of it.
@named F = Se(sin(t))
@named m = Inertia(I = 1.0)
@named k = Compliance(C = 1.0)
@named b = Resistance(R = 0.5)
@named model = Junction1(F => :pout, m, k, b)

# The equations in terms of the states, the spring's displacement `q` and the mass's momentum `p`. Namespaced names such as `k₊q` show which element each variable belongs to:
bg_equations(model)

# ## Renaming
#
# Map any variable or parameter to the name to show. Renaming a variable renames its derivative too, unless the derivative has its own entry, here `ẋ` for the velocity of the spring's end:
@variables x(t) ẋ(t) p(t)
@parameters M c C
eqs = bg_equations(model; rename = Dict(model.k.q => x, D(model.k.q) => ẋ, model.m.p => p,
                                        model.m.I => M, model.b.R => c, model.k.C => C))

# ## LaTeX
latexify(eqs)

# ## Renamed variables in a solution
#
# `bg_equations` only changes how the equations are shown. To read a renamed variable from a solution, add an alias equation to the model; MTK keeps it as an observed variable:
@named aliased = System([x ~ model.k.q], t; systems = [model])
sol = solve(ODEProblem(mtkcompile(aliased), [], (0.0, 10.0)), Tsit5())
sol[x][end]

# A solution can be indexed by one variable at a time. To plot or compute a product such as the power `e⋅f` on a bond, sample each factor and multiply:
power = sol[model.b.e] .* sol[model.b.f]
power[end]
