# # [MSD system with N-DoF](@id msd-ndof)
# 
# One of the advantages of using a scripting approach to model dynamic systems is the ability to create and simulate the system algorithmically. To showcase this capability using the Bond Graph Toolkit, let us consider an example of a mass-spring-damper (MSD) system with $n$ degrees of freedom. We will use a `for` loop to vary the number of masses. The n-degree-of-freedom (n-DoF) MSD system is depicted below:
# 
# TODO: Update the image
# ![Mass-Spring-Damper System](../assets/block_2dof.png)
# 
# The following differential equation governs the motion of this n-DoF harmonic oscillator:
# 
# TODO: Update the equation
# $m\,\frac{d^2 x}{dt^2} + c\,\frac{d^2 x}{dt^2} + k\,x = F(t)$,
# 
# where $x_n$ is the position, $m_n$ is the mass, $c_n$ is the damping coefficient, $k_n$ is the spring stiffness, and $F(t)$ is an external force applied to the system.
# 
# To model the n-DoF system using the bond graph method, we need to connect inertance (I), compliance (C), and resistance (R) elements to different 1-junctions and 0-junctions, as shown in the bond graph diagram:
# 
# TODO: Generate the bond graph for the n-DoF system
# ![Bond Graph Representation](../assets/bg_2dof.png)
# 
# In this representation, the inertance element (I) corresponds to the mass (m) in the mechanical system, the compliance element (C) corresponds to the spring constant (k), and the resistance element (R) represents the damping coefficient (c). For simplicity, we assume they are equal for both degrees of freedom.
# 
# ## Bond Graph Toolkit
# 
# We need to import the Bond Graph Toolkit module and the `DifferentialEquations.jl` package to model the system and solve the resultant ODE. Additionally, we will import the independent variable `t` from the Bond Graph Toolkit for defining a custom forcing term.

using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# ## Building the model
# 
# ### Setting the elements
# 
# Before defining a system, we need to define the single port elements of the system. Then, we define the mass, spring, and damper elements as follows:
@named m = Mass(m = 1.0)      # (kg)
@named s = Spring(k = 1.0)    # (N/m)
@named d = Damper(c = 1.0)    # (N⋅m⋅s)

# ### Setting the junctions and subsystems
# 
# Similarly to the [2-DoF](@ref msd-2dof) example, in the n-DoF system, we need to define multiple 1-junctions and 0-junctions. To simplify the process, we can define them separately and later connect them accordingly. In the bond graph model of the n-DoF system, there is a section that repeats across the model, a spring, damper, and one-junction. To avoid redundancy, we can create this section as a subsystem and reuse it in other parts of the system where it is needed.

# Firstly we define the ground-connected mass ($m_1$) as a subsystem.
@named m1 = Junction1([-1, m], [-1, d], [-1, s])

# Then, we create the n-DoF spring and damper subsystem.
@named sd = Junction1([-1, s], [-1, d])

# Set the number of degrees of freedom. For this case, it should be greater or equal to $3$.
n = 4
# Then, we build the model using a `for` loop. The subsystems are stored inside the `Matrix` of subsystems (`M`), while the connections between the subsystems are stored inside the `Vector` of connections (`cons`).
M, cons = [], []
for i = 2:n
    ## Generate the subsystems of one DoF
    j0 = Junction0(sd; name = Symbol("m" * string(i) * "_j0"))
    j1 = Junction1([-1, m]; name = Symbol("m" * string(i) * "_j1"))
    push!(M, [j0, j1])
    ## Connect the subsystems of one DoF
    push!(cons, connect(j0.power, j1.power))
end

# ### Connect the subsystems

# We need now to iterate over the subsystem connecting them.
push!(cons, [connect(M[i][2].power, M[i+1][1].power) for i = 1:(length(M)-1)]...)

# Finally, we need to connect the ground-connected mass ($m_1$) to the second subsystem.
push!(cons, connect(m1.power, M[1][1].power))

# ### Build the system 

# We flatten the subsystem `Matrix` into a `Vector` for simplicity.
M_flat = vcat(M...);

# We build the system by creating an `ODESystem` with the connections and setting the independent variable as t.
@named mdl = ODESystem(cons, t)
# Then, we need to add the subsystems to the model.
mdl = compose(mdl, m1, M_flat...)

# ## Analysing the model
# We can visualize the bond-graph connections by generating a graph plot of the model.
generate_graph(mdl)

# We can obtain the bond-graph equations using the code below.
equations(expand_connections(mdl))

# The model's equations above were not simplified and represented the bond-graphs equations directly derived from the connections. Then, we simplify the DAE above and obtain the ODE. The simplified set of equations is required for using `solve` from `DifferentialEquations.jl`.
@named sys = simplifysys(mdl);

# Print system states.
states(sys)

# Print system parameters.
parameters(sys)

# Print the system's simplified equations.
equations(sys)

# ## Simulate the system
# 
# ### Unforced case
# In this case, we will consider an initial displacement of $1.0 (m)$ on the second mass of the system.

# Define the simulation time.
tspan = (0, 10)

# When defining the ODEProblem as in the ModelingToolkit, we can define the initial value of the states.
prob = ODEProblem(sys, [M[1][2].m.power.f => 0.0, M[1][1].sd.s.q => 1.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel="Time", ylabel="Amplitude")