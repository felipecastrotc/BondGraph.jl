# # [Harmonic Oscillator](@id msd-1dof)
# 
# ## Introduction
# 
# The harmonic oscillator exemplifies dynamic systems with oscillatory behaviors. When considering the mechanical domain, this system manifests as a mass-spring-damper system, as illustrated below:
# 
# ```@raw html
# <picture style="display: flex; justify-content: center; align-items: center;">
#   <source srcset="../../assets/tikz/fbd_1dof-dark.svg" media="(prefers-color-scheme: dark)">
#   <img src="../../assets/tikz/fbd_1dof.svg">
# </picture>
# ```
# 
# The governing differential equation for this oscillator is:
# 
# $ m\,\frac{d^2 x}{dt^2} + c\,\frac{dx}{dt} + k\,x = F(t), $
# 
# Where:
# - $ x $ denotes the position,
# - $ m $ is the mass,
# - $ c $ represents the damping coefficient,
# - $ k $ is the spring stiffness, and
# - $ F(t) $ is the external force applied.
# 
# By utilizing the bond graph approach, this system can be modeled by connecting inertance (I), compliance (C), and resistance (R) elements to a one-junction:

# ```@raw html
# <picture style="display: flex; justify-content: center; align-items: center;">
#   <source srcset="../../assets/tikz/bg_1dof-dark.svg" media="(prefers-color-scheme: dark)">
#   <img src="../../assets/tikz/bg_1dof.svg">
# </picture>
# ```

# In this representation, the inertance element (I) corresponds to the mass (m) in the mechanical system, the compliance element (C) corresponds to the spring constant (k), and the resistance element (R) represents the damping coefficient (c).

# ## Preliminary Steps

# Firstly, we need to import the bond-graph toolkit module to model the system and the `DifferentialEquations.jl` to solve the resultant ODE. Also, we import the independent variable `t` from the bond-graph toolkit to define a custom forcing term.

using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# ## System Element Definitions
# 
# Before defining a system, we need to define the single port elements of the system. Then, we define the mass, spring, and damper elements as follows:
@named m = Mass(m = 1.0);      # (kg)
@named s = Spring(k = 1.0);    # (N/m)
@named d = Damper(c = 1.0);    # (N⋅m⋅s)

# The forcing term can be set as an expression or even a Julia function. In this example, we will first initialize a generic forcing term that will be defined before solving the ODE.
@parameters F(t)
@named f = Se(F)

# ## Model Construction
# 
# Then, we build the system by passing the one-port elements to the 1-junction function as arguments.
@named mdl = Junction1(m, d, s, [-1, f]);

# ## Analysing the model
# We can visualize the bond-graph connections by generating a graph plot of the model.
generate_graph(mdl)

# We can obtain the bond-graph equations using the code below.
equations(expand_connections(mdl))

# The model's equations above were not simplified and represent the bond-graphs equations directly derived from the connections. Then, we simplify the DAE above and obtain the ODE. The simplified set of equations is required for using `solve` from `DifferentialEquations.jl`.
@named sys = simplifysys(mdl)

# Print the system's states.
states(sys)

# Print the system's parameters.
parameters(sys)

# Print the system's simplified equations.
equations(sys)

# By comparison, we can see that the equations above are the same as the canonical form of the harmonic oscillator presented at the beginning.

# We can generate LaTeX code from the equations.
latexify(equations(sys))

# ## System Simulations
# 
# ### Unforced case
# Now, we will only the unforced case. Therefore, $F(t) = 0$. Thus, we update the force function on the system.
sys_unforced = substitute(sys, Dict(f.F => 0.0))

# Define the simulation time.
tspan = (0, 10)

# When defining the ODEProblem as in the ModelingToolkit, we can define the initial value of the states.
prob = ODEProblem(sys_unforced, [m.power.f => 0.0, s.q => 1.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")

# ### Forced case
# Now, we will consider the forced case, where $F(t) = sin(t)$. Thus, we update the force function on the system.
sys_forced = substitute(sys, Dict(f.F => sin(t)))

# Define the simulation time.
tspan = (0, 10)

# When defining the ODEProblem, as in the ModelingToolkit, we can define the initial value of the states.
prob = ODEProblem(sys_forced, [m.power.f => 0.0, s.q => 0.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")
