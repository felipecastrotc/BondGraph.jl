# # [Quick Example](@id rc-circuit)
#
# In this quick example, we'll showcase how to model a basic RC circuit using the bond graph method. An RC circuit comprises a resistor (R) and a capacitor (C). Below is the visual representation of such a circuit:
#
# ![RC Circuit](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a4/Discharging_capacitor.svg/400px-Discharging_capacitor.svg.png)
#
# The capacitor discharges through the resistor. Both share the voltage $V = q/C$, and Kirchhoff's current law says the current leaving the capacitor, $-\dot{q}$, is the current through the resistor, $V/R$. So the charge $q$ on the capacitor obeys
#
# $R \frac{dq}{dt} + \frac{q}{C} = 0, \qquad q(t) = q_0\, e^{-t/RC}.$
#
# ### Bond Graph Representation
#
# The resistor and the capacitor share one voltage, so in the bond graph they hang off a 0-junction. As usual for passive elements, both bonds point from the junction into the element (Karnopp, Margolis & Rosenberg, §3.1):
#
# ```@raw html
# <div class="docs-light-only" style="text-align: center;">
# ```
# ![Bond graph](../assets/tikz/bg_rc.svg)
# ```@raw html
# </div><div class="docs-dark-only" style="text-align: center;">
# ```
# ![Bond graph](../assets/tikz/bg_rc-dark.svg)
# ```@raw html
# </div>
# ```
#
# ### Using the Bond Graph Toolkit
#
# Follow the steps below to model, obtain equations, and simulate the RC circuit using the Bond Graph Toolkit:

## Import the required libraries
using BondGraphToolkit
using OrdinaryDiffEqTsit5
using Plots

## Define the elements of the RC circuit
@named C = Compliance(C = 1.0, q = 1.0)   # capacitance (F), initial charge (C)
@named R = Resistance(R = 1.0)            # resistance (Ohm)

## Bond both elements to a 0-junction; a bare element `x` means the bond points out of the junction, into x
@named model = Junction0(R, C)
## Compile the model into its simplest form
sys = mtkcompile(model)

## Define and solve the ODE problem over 10 s
prob = ODEProblem(sys, [], (0.0, 10.0))
sol = solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-10)
## Plot the charge on the capacitor against the hand-derived q(t) = q₀ exp(-t/RC), with q₀ = R = C = 1
plot(sol, idxs = [model.C.q], label = "Bond graph", xlabel = "Time", ylabel = "Charge")
plot!(τ -> exp(-τ), 0, 10, label = "Hand-derived", linestyle = :dash)

# The largest difference between the two on a common time grid:
tgrid = 0:0.01:10
maximum(abs, sol(tgrid; idxs = model.C.q).u - exp.(-tgrid))

# To inspect the system, the following functions can be employed:
#
# - **State Inspection:**
unknowns(sys)

# - **Parameter Inspection:**
parameters(sys)

# - **Equation Inspection:**
equations(sys)

# ### Code Explanation
#
# 1. **Elements.** `Compliance` and `Resistance` are the bond-graph C and R elements. Their laws are `q̇ = f`, `e = q / C` and `e = R f`, written in the variables of the element's bond.
#
# 2. **Junction.** `Junction0` makes every bond share one effort (the voltage) and sums the flows with the sign of each bond's half-arrow: `+` into the junction, `-` out of it. Writing `R` draws the bond from the junction into `R`; `x => :pout` would draw it from `x` into the junction.
#
# 3. **Simulation.** A bond-graph model is an ordinary ModelingToolkit system, so `mtkcompile`, `ODEProblem` and `solve` work as for any MTK model. The initial charge comes from the `q = 1.0` keyword of `Compliance`.
