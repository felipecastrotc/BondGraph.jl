# # [Mass-Spring-Damper system 2-DoF](@id msd-2dof)
# 
# ## Introduction
# 
# A more challenging yet simple dynamic system than the [harmonic oscillator](@ref msd-1dof) is the two degrees of freedom mass-spring-damper system, illustrated below: 
# 
# ```@raw html
# <picture style="display: flex; justify-content: center; align-items: center;">
#   <source srcset="../../assets/tikz/fbd_2dof-dark.svg" media="(prefers-color-scheme: dark)">
#   <img src="../../assets/tikz/fbd_2dof.svg">
# </picture>
# ```
# 
# The governing differential equation for the two-degree of freedom harmonic oscillator is:
# 
# TODO: update the equation
# $m\,\frac{d^2 x}{dt^2} + c\,\frac{d^2 x}{dt^2} + k\,x = F(t)$,
# where:
# - $ x $ is the position,
# - $ m $ is the mass,
# - $ c $ is the damping coefficient,
# - $ k $ is the spring stiffness,
# - $ F(t) $ is an external force applied to the system.
# 
# To model the 2-DoF system using the bond graph method, we need to connect the inertance (I), compliance (C), and resistance (R) elements to different 1-junctions and 0-junctions, as depicted in the bond graph model below:
# 
# ```@raw html
# <picture style="display: flex; justify-content: center; align-items: center;">
#   <source srcset="../../assets/tikz/bg_2dof-dark.svg" media="(prefers-color-scheme: dark)">
#   <img src="../../assets/tikz/bg_2dof.svg">
# </picture>
# ```
# In the bond graph:
# - The inertance (I) element corresponds to mass \( m \),
# - The compliance (C) element corresponds to spring constant \( k \),
# - The resistance (R) element represents the damping coefficient \( c \).

# For simplicity, we assume equal values for both degrees of freedom.

# ## Preliminary Steps
# 
# Firstly, we need to import the bond-graph toolkit module to model the system and the `DifferentialEquations.jl` to solve the resultant ODE. Also, we import the independent variable `t` from the bond-graph toolkit for defining a custom forcing term.

using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# ## System Element Definitions
# 
# Before defining a system, we need to define the single port elements of the system. Then, we define the mass, spring, and damper elements as follows:
@named m = Mass(m=1.0);      # (kg)
@named s = Spring(k=1.0);    # (N/m)
@named d = Damper(c=1.0);    # (N⋅m⋅s)

# The forcing term can be set as an expression or even a Julia function. In this example, we will first intialize a generic forcing term that will be defined before solving the ODE.
@parameters F(t)
@named f = Se(F)

# ### Setting the junctions and subsystems
# 
# Differently from the [1-DoF](@ref msd-1dof), where we could define the entire system with only a 1-junction. For the 2-DoF system we are required to define multiple 1-junctions and one 0-junction. Thus, we define them separately and later we will connect them accordingly.
# 
# Firstly we define the ground connected mass ($m_1$).
@named m1_j1 = Junction1([-1, m], [-1, d], [-1, s])

# Then, we define the $2^{nd}$ DoF mass (($m_2$)) with two 1-junctions and one 0-junction.
@named m2_j1 = Junction1([-1, m], f)
@named m2_j0 = Junction0()
@named m2_sd = Junction1([-1, s], [-1, d])

# ## Model Construction
# 
# ### Connect the subsystems

# Define the connections between the junctions.
cons = [
    connect(m2_j0.power, m2_j1.power),
    connect(m2_j0.power, m2_sd.power),
    connect(m1_j1.power, m2_j0.power),
]

# ### Build the system 
# 
# We build the system by first creating a `ODESystem` with junctions connections and setting the independent variable as t.
@named mdl = ODESystem(cons, t)
# Then, we need to add the junctions subsystems to the model
mdl = compose(mdl, m1_j1, m2_j0, m2_j1, m2_sd)

# ## Analysing the model
# 
#  We can visualize the bond-graph connections by generating a graph plot of the model.
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

# We can generate LaTeX code from the equations
latexify(equations(sys))

# ## System Simulations
# 
# ### Unforced case
# Now, we will only the unforced case. Therefore, $F(t) = 0$. Thus, we update the force function on the system.
sys_unforced = substitute(sys, Dict(m2_j1.f.F => 0.0))

# Define the simulation time.
tspan = (0, 10)

# When defining the ODEProblem as in the ModelingToolkit, we can define the initial value of the states.
prob = ODEProblem(sys_unforced, [m2_j1.m.power.f => 0.0, m2_sd.s.q => 1.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel="Time", ylabel="Amplitude")

# ### Forced case
# Now, we will consider the forced case, where $F(t) = sin(t)$. Thus, we update the force function on the system.
sys_forced = substitute(sys, Dict(m2_j1.f.F => sin(t)))

# Define the simulation time.
tspan = (0, 10)

# When defining the ODEProblem, as in the ModelingToolkit, we can define the initial value of the states.
prob = ODEProblem(sys_forced, [m2_j1.m.power.f => 0.0, m2_sd.s.q => 0.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel="Time", ylabel="Amplitude")