# # [Harmonic Oscillator](@id msd-1dof)
# 
# The harmonic oscillator is a classic example of a dynamic system characterized by oscillatory behavior. In the mechanical domain, it can be represented as a mass-spring-damper system, as illustrated below:
# 
# ![Mass-Spring-Damper System](../assets/block_1dof.png)
# 
# The following differential equation governs the motion of the harmonic oscillator:
# 
# $m\,\frac{d^2 x}{dt^2} + c\,\frac{d^2 x}{dt^2} + k\,x = F(t), $
# 
# where $x$ is the position, $m$ is the mass, $c$ is the damping coefficient, $k$ is the spring stiffness, and $F(t)$ is an external force applied to the system.
# 
# To model the harmonic oscillator using the bond graph method, we can connect the inertance (I), compliance (C), and resistance (R) elements to a 1-junction, as depicted in the bond graph diagram:

# ![Bond Graph Representation](../assets/bg_1dof.png)

# In this representation, the inertance element (I) corresponds to the mass (m) in the mechanical system, the compliance element (C) corresponds to the spring constant (k), and the resistance element (R) represents the damping coefficient (c).

# ## Bond-graph toolkit

# Firstly, we need to import the bond-graph toolkit module to model the system and the `DifferentialEquations.jl` to solve the resultant ODE. Also, we import the independent variable `t` from the bond-graph toolkit to define a custom forcing term.

using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# ## Building the model
# 
# Before defining a system, we need to define the single port elements of the system. Then, we define the mass, spring, and damper elements as follows:
@named m = Mass(m = 1.0);      # (kg)
@named s = Spring(k = 1.0);    # (N/m)
@named d = Damper(c = 1.0);    # (N⋅m⋅s)

# The forcing term can be set as an expression or even a Julia function. In this example, we will first intialize a generic forcing term that will be defined before solving the ODE.
@parameters F(t)
@named f = Se(F)

# Then, we build the system by simply passing the one-port elements to the 1-junction function as arguments.
@named mdl = Junction1(m, d, s, [-1, f]);

# ## Analysing the model
# We can visualize the bond-graph connections by generating a graph plot of the model.
generate_graph(mdl)

# We can obtain the bond-graph equations using the code below.
equations(expand_connections(mdl))

# The equations of the model above were not simplified and represent the bond-graphs equations directly derived from the connections. Then, we simplify the DAE above and obtain the ODE. The simplified set of equations is required for using `solve` from `DifferentialEquations.jl`
@named sys = simplifysys(mdl)

# Print system states
states(sys)

# Print system parameters
parameters(sys)

# Print the system simplified equations
equations(sys)

# By comparison, we can see that the equations above are the same as the canonical form of the harmonic oscillator presented at the beginning.

# We can generate LaTeX code from the equations
latexify(equations(sys))

# ## Simulate the system
# 
# ### Unforced case
# Now, we are going to consider only the unforced case. Therefore, $F(t) = 0$. Thus, we update the force function on the system.
sys_unforced = substitute(sys, Dict(f.F => 0.0))

# Define the simulation time
tspan = (0, 10)

# When defining the ODEProblem as in the ModelingToolkit, we can define the initial value of the states.
prob = ODEProblem(sys_unforced, [m.power.f => 0.0, s.q => 1.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")

# ### Forced case
# Now, we are going to consider the forced case, where $F(t) = sin(t)$. Thus, we update the force function on the system.
sys_forced = substitute(sys, Dict(f.F => sin(t)))

# Define the simulation time
tspan = (0, 10)

# When defining the ODEProblem as in the ModelingToolkit we can define the initial value of the states.
prob = ODEProblem(sys_forced, [m.power.f => 0.0, s.q => 0.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")
