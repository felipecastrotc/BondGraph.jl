# # [Quick Example](@id rc-circuit)
#
# In this quick example, we'll showcase how to model a basic RC circuit using the bond graph method. An RC circuit comprises a resistor (R) and a capacitor (C). Below is the visual representation of such a circuit:
#
# ![RC Circuit](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a4/Discharging_capacitor.svg/400px-Discharging_capacitor.svg.png)
#
# From Kirchhoff's law, the linear differential equation representing the charge on the capacitor over time is given by:
#
# $C \frac{dV}{dt} + \frac{V}{R} = 0$
#
# ### Bond Graph Representation
#
# The resistor and the capacitor share one voltage, so in the bond graph they hang off a 0-junction. As usual for passive elements, both bonds point from the junction into the element (Karnopp, Margolis & Rosenberg, §3.1):
#
# ```@raw html
# <picture style="display: flex; justify-content: center; align-items: center;">
#   <source srcset="../../assets/tikz/bg_rc-dark.svg" media="(prefers-color-scheme: dark)">
#   <img src="../../assets/tikz/bg_rc.svg">
# </picture>
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

## Bond both elements to a 0-junction; `[-1, x]` means the bond points out of the junction, into x
@named model = Junction0([-1, R], [-1, C])
## Compile the model into its simplest form
sys = mtkcompile(model)

## Define and solve the ODE problem over 10 s
prob = ODEProblem(sys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
## Plot the charge on the capacitor
plot(sol, idxs = [model.C.q], xlabel = "Time", ylabel = "Charge")

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
# 2. **Junction.** `Junction0` makes every bond share one effort (the voltage) and sums the flows with the sign of each bond's half-arrow: `+` into the junction, `-` out of it. Writing `[-1, R]` draws the bond from the junction into `R`.
#
# 3. **Simulation.** A bond-graph model is an ordinary ModelingToolkit system, so `mtkcompile`, `ODEProblem` and `solve` work as for any MTK model. The initial charge comes from the `q = 1.0` keyword of `Compliance`.
