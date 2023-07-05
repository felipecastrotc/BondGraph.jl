# # [Duffing Equation](@id msd-duffing)
# 
# The Duffing equation is a nonlinear second-order ordinary differential equation commonly used to model various physical systems exhibiting nonlinear behavior, such as electrical circuits, mechanical systems, and biological systems. From the mechanical domain point of view, it is an extension of the simple mass-spring-damper system that incorporates a cubic spring. Despite its simple form, the Duffing equation exhibits chaotic behavior. The equation is given by:
# 
# $m\,\frac{d^2 x}{dt^2} + c\,\frac{d x}{dt} + k\,x + k_3\,x^3 = F(t),$
# 
# where $x$ represents the position, $m$ is the mass, $c$ is the damping coefficient, $k$ is the spring stiffness, $k_3$ is the cubic spring stiffness, and $F(t)$ denotes an external force applied to the system.
# 
# In standard bond graphs, the capacitance element is linear, while for the Duffing equation, we need a nonlinear capacitance for the cubic spring. Therefore, we need to introduce a nonlinear capacitance element to model the Duffing equation using the Bond Graph Toolkit. The library allows users to define custom elements and utilize the [`Power`](@ref) connector as interface.
# 
# # Bond Graph Toolkit
# 
# First, we import the Bond Graph Toolkit module to model the system and the `DifferentialEquations.jl` package to solve the resulting ordinary differential equation (ODE). Additionally, we import the independent variable `t` from the Bond Graph Toolkit, which will be used to define a custom forcing term.

using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# ## Building the model
# 
# ### Custom capacitance
# 
# We now define a custom element called `Spring3` to include the non-linearity of the Duffing equation. The element is defined by the function below.

function Spring3(; name, k = 1.0, x = 0.0)
    ## Initialize the Power connector
    @named power = Power()

    ## Initialize the displacement state with its initial value `x`
    @variables q(t) = x
    ## Set the equation parameters
    ps = @parameters C = 1 / k

    ## Define Spring3 element equations with the non-lineariyy of the cubic spring
    eqs = [
        power.e ~ q^3 / C
        D(q) ~ power.f
    ]

    ## Build the system with the equations, state and connector.
    compose(ODESystem(eqs, t, [q], ps; name = name), power)
end

# The `Spring3` function creates a bond graph subsystem representing the cubic spring. You can specify the following parameters:
# 
# - `name`: The name of the `Spring3` element (default: "").
# - `k`: The cubic spring stiffness constant (default: 1.0).
# - `x`: The initial displacement of the spring (default: 0.0).
# 
# ### Setting the elements
# 
# Now we initialize and define the cubic spring custom element and the remaining single port elements of the system (mass, spring, and damper) as follows:

# - The parameters values
α = 1.0
β = 5.0
δ = 0.02
γ = 8.0
ω = 0.5

# - Define the elements
@named m = Mass(m = 1.0)
@named s = Spring(k = α)
@named s3 = Spring3(k = β)
@named d = Damper(c = δ)

# The forcing term can be set as an expression or even a Julia function. In this example, we will first intialize a generic forcing term that will be defined before solving the ODE.
@parameters F(t)
@named f = Se(F)

# Then, we build the system by simply passing the one-port elements to the 1-junction function as arguments.
@named mdl = Junction1(m, d, s, s3, [-1, f])

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
prob = ODEProblem(sys_unforced, [m.power.f => 0.0, s.q => 1.0, s3.q => 1.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")

# ### Forced case
# Now, we are going to consider the forced case, where $F(t) = sin(t)$. Thus, we update the force function on the system.
sys_forced = substitute(sys, Dict(f.F => γ * cos(ω * t)))

# Define the simulation time
tspan = (0, 10)

# When defining the ODEProblem as in the ModelingToolkit we can define the initial value of the states.
prob = ODEProblem(sys_forced, [m.power.f => 0.0, s.q => 0.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")