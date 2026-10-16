# # [Transformers and Gyrators](@id two-ports)
#
# A transformer (TF) and a gyrator (GY) are two-port elements that conserve power (Karnopp, Margolis & Rosenberg, §3.2). Their port `p1` is the bond into the element and `p2` the bond out of it:
#
# - TF with modulus `m`: `e₁ = m e₂`, `f₂ = m f₁` (a lever, a gear pair, a hydraulic ram);
# - GY with modulus `r`: `e₁ = r f₂`, `e₂ = r f₁` (a DC motor, a voice coil).
#
# Here a constant force drives a mass through a lever, first with a fixed ratio, then with a ratio that changes in time, and finally with a transformer written by hand.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t
using OrdinaryDiffEqTsit5
using Plots

# ## A lever: constant transformer
#
# The force source acts on a 1-junction, whose free port `out` bonds into the lever's port 1. The lever's port 2 bonds straight into the mass: a plain `connect` between a two-port and an element is a bond pointing into the element.
function lever(tf)
    @named F = Se(3.0)
    @named mass = Inertia(I = 1.0)
    @named src = Junction1(F, [-1, :out])
    System([connect(src.out, tf.p1), connect(tf.p2, mass.power)], t;
           name = :model, systems = [src, tf, mass])
end

@named tf = TF(m = 2.0)
model = lever(tf)
sol = solve(ODEProblem(mtkcompile(model), [], (0.0, 5.0)), Tsit5())

# The mass feels `e₂ = e₁ / m = 1.5`, so its momentum grows as `p = 1.5 t`.
plot(sol, idxs = [model.mass.p], xlabel = "Time", ylabel = "Momentum")

# ## Modulated transformer
#
# Passing an expression as the modulus gives a modulated transformer (MTF). Here the lever arm grows in time, `m = 1 + t`, so the force on the mass decays and `p = 3 ln(1 + t)`.
@named tf = TF(m = 1 + t)
model = lever(tf)
sol = solve(ODEProblem(mtkcompile(model), [], (0.0, 5.0)), Tsit5())
plot(sol, idxs = [model.mass.p], xlabel = "Time", ylabel = "Momentum")

# The modulus can also read the model's own variables, as in the [DC motor](@ref dc-motor) example. Build such a two-port with `@named`: it scopes the other subsystems' variables to the parent model.
#
# ## A custom two-port
#
# `@twoport` defines a new two-port element. It adds the ports `p1` (bond in) and `p2` (bond out) and gives you the bonds' own variables `e1, f1, e2, f2`, each flow positive along its half-arrow, so the law is written as in the book. This transformer is the same as the built-in `TF`:
@twoport function MyTF(; name, m)
    @parameters m = m
    [e1 ~ m * e2, f2 ~ m * f1]
end

@named tf = MyTF(m = 2.0)
model = lever(tf)
sol = solve(ODEProblem(mtkcompile(model), [], (0.0, 5.0)), Tsit5())
plot(sol, idxs = [model.mass.p], xlabel = "Time", ylabel = "Momentum")
