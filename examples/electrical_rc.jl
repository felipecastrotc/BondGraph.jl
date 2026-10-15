# ## [Quick Example](@id rc-circuit)
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
# To model the RC circuit using the bond graph method, we connect the compliance (C) and resistance (R) elements to a zero-junction, as shown in the bond graph diagram:
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
using BondGraph
using DifferentialEquations
using Plots

## Define the elements of the RC circuit
@named C = Spring(k=1 / 1.0)   # Electric capacitance (F)
@named R = Damper(c=1.0)       # Electric resistance (Ohm)

## Create the RC circuit model by connecting the elements
@named model = Junction0(R, C)
## Simplify the system equations
@named sys = simplifysys(model)

## Set the time span for simulation
tspan = (0, 10)
## Define the initial conditions
initial_conditions = [R.power.f => 0.0, C.q => 1.0]

## Define the ODE problem
prob = ODEProblem(sys, initial_conditions, tspan)
## Solve the ODE problem
sol = solve(prob)
## Plot the results
plot(sol, xlabel="Time", ylabel="Amplitude")

# To visualize the bond graph diagram of the model, utilize the `generate_graph` function:

generate_graph(model)

# To inspect the system, the following functions can be employed:
# 
# - **State Inspection:** 
states(sys)

# - **Parameter Inspection:** 
parameters(sys)

# - **Equation Inspection:** 
equations(sys)

# ### Code Explanation
# 
# For clarity, let's break down the code snippets step by step:
# 
# 1. **Library Importation:** 
# 
#    The libraries `BondGraph`, `DifferentialEquations`, and `Plots` are imported. These are essential for bond graph modeling, differential equation solving, and data visualization respectively.
# 
#    ```julia
#    using BondGraph
#    using DifferentialEquations
#    using Plots
#    ```
# 
# 2. **Defining Circuit Elements:**
# 
#    The elements of the RC circuit, specifically the compliance `C` and the resistance `R`, are defined. For this, we employ the `Spring` and `Damper` functions from the `BondGraphToolkit`.
# 
#    ```julia
#    @named C = Spring(k=1 / 1.0)  
#    @named R = Damper(c=1.0)      
#    ```
# 
# 3. **Model Creation:**
# 
#    The RC circuit model is constructed by connecting the `R` and `C` elements to a zero-junction using `Junction0`. To streamline the model, the `simplifysys` function is applied, and the result is assigned to the `sys` variable.
# 
#    ```julia
#    @named model = Junction0(R, C)
#    @named sys = simplifysys(model)
#    ```
# 
# 4. **Simulation Setup:**
# 
#    The time span for the simulation is defined as `(0, 10)`. The initial conditions dictate that the capacitor `C` starts with a charge of $1.0$, while the flow through resistor `R` is initialized at $0.0$.
# 
#    ```julia
#    tspan = (0, 10)
#    initial_conditions = [R.power.f => 0.0, C.q => 1.0]
#    ```
# 
# 5. **Solving the Differential Equations:**
# 
#    The differential equations that describe the RC circuit are set up using the `ODEProblem` function. The function's arguments include the system `sys`, the initial conditions, and the simulation time span. The problem is then solved using the `solve` function.
# 
#    ```julia
#    prob = ODEProblem(sys, initial_conditions, tspan)
#    sol = solve(prob)
#    ```
# 
# 6. **Visualization:**
# 
#    The simulation results are depicted using the `plot` function. The graph's x-axis represents time, while the y-axis indicates the amplitude of the system variables.
# 
#    ```julia
#    plot(sol, xlabel="Time", ylabel="Amplitude")
#    ```