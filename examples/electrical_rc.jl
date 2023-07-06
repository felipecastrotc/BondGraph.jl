# ## [Quick Example](@id rc-circuit)
#
# Let's explore the capabilities of the library through a quick example of modeling a basic dynamic system, an RC circuit. The RC circuit is a fundamental electrical circuit consisting of a resistor (R) and a capacitor (C), as depicted in the figure below:
#
# ![RC Circuit](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a4/Discharging_capacitor.svg/400px-Discharging_capacitor.svg.png)
#
# By using the Kirchhoff's law we have the following linear differential equation that describes the charge on the capacitor along the time.
#
# $C \frac{dV}{dt} + \frac{V}{R} = 0$
#
# To model the RC circuit using the bond graph method, we can connect the capacitance (C) and resistance (R) elements to a zero-junction, as shown in the bond graph diagram:
#
# ![Bond Graph Representation](../assets/bg_1dof.png)
#
# ### Bond Graph Toolkit
#
# We can use the Bond Graph Toolkit to model and get the RC circuit equations and then solve the system, as presented in the code below.

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

# Moreover, we can generate the bond graph diagram of the model using the `generate_graph` function. 

generate_graph(model)

# Print the system states
states(sys)

# Print the system parameters
parameters(sys)

# Print the simplified equations of the system
equations(sys)

# In the code snippets, we start by importing the libraries `BondGraph`, `DifferentialEquations`, and `Plots`. Then, ee define the elements of the RC circuit, namely the capacitance `C` and the resistance `R`, using the functions `Spring` and `Damper` from the `BondGraphToolkit` package.

# Next, we create the RC circuit model by connecting the elements `R` and `C` to a zero-junction using the `Junction0` function. We use the `simplifysys` function to obtain a more concise representation of the system, which is stored in the variable `sys`.

# To simulate the system, We first define the time span for the simulation using the tuple `(0, 10)`. Then, we specify the initial conditions of the system in the `initial_conditions` dictionary, where we set the initial charge on the capacitor `C` to $1.0$ and the initial flow through the resistor `R` to $0.0$.

# Using the `ODEProblem` function from `DifferentialEquations`, we define the ODE problem to be solved by passing the system `sys`, the initial conditions, and the time span. We solve the ODE problem using the `solve` function.

# Finally, we visualize the results by plotting the solution (`sol`) using the `plot` function from `Plots`. The x-axis represents time, and the y-axis represents the amplitude of the variables in the system.
