# # [Duffing Equation](@id msd-duffing)
# 
# ## Introduction
# 
# The Duffing equation represents a nonlinear second-order ordinary differential equation. Widely employed in the modeling of physical systems exhibiting nonlinear behavior, such as electrical circuits, mechanical setups, and even biological systems, this equation stands out for its rich dynamics. In the mechanical domain, the Duffing equation is an enhanced version of the conventional mass-spring-damper system, incorporating a cubic spring. The equation is mathematically represented as:
# 
# $ m\,\frac{d^2 x}{dt^2} + c\,\frac{d x}{dt} + k\,x + k_3\,x^3 = F(t), $
# 
# where:
# 
# - $ x $: Position
# - $ m $: Mass
# - $ c $: Damping coefficient
# - $ k $: Spring stiffness
# - $ k_3 $: Cubic spring stiffness
# - $ F(t) $: External force
# 
# For visualization, the system's free-body diagram is depicted below:
# 
# ```@raw html
# <picture style="display: flex; justify-content: center; align-items: center;">
#   <source srcset="../../assets/tikz/fbd_duffing-dark.svg" media="(prefers-color-scheme: dark)">
#   <img src="../../assets/tikz/fbd_duffing.svg">
# </picture>
# ```
# 
# Its corresponding representation in bond graph notation is:
# 
# ```@raw html
# <picture style="display: flex; justify-content: center; align-items: center;">
#   <source srcset="../../assets/tikz/bg_duffing-dark.svg" media="(prefers-color-scheme: dark)">
#   <img src="../../assets/tikz/bg_duffing.svg">
# </picture>
# ```
# 
# Given the Duffing equation's nonlinear nature, standard bond graphs require a modification. Specifically, a nonlinear compliance element is introduced to represent the cubic spring, all achieved through the Bond Graph Toolkit, which offers capabilities to customize and incorporate such elements through the `Power` connector.
# 
# ## Preliminary Steps
# 
# First, we import the Bond Graph Toolkit module to model the system and the `DifferentialEquations.jl` package to solve the resulting ordinary differential equation (ODE). Additionally, we import the independent variable `t` from the Bond Graph Toolkit, which will be used to define a custom forcing term.

using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# ## System Element Definitions
# 
# ### Custom Compliance for Cubic Spring
# 
# To account for the Duffing equation's non-linearity we define a custom element called `Spring3`. The custom element definition can be simplified using the macro `@oneport` as shown below.

@oneport function Spring3(; name, k = 1.0, x = 0.0)

    ## Initialize the displacement state with its initial value `x`
    @variables q(t) = x
    ## Set the equation parameters
    @parameters C = 1 / k

    ## Define Spring3 element equations with the non-linearity of the cubic spring
    [
        power.e ~ q^3 / C
        D(q) ~ power.f
    ]
end

# The `Spring3` function creates a bond graph subsystem representing the cubic spring, which can be customized via:
# 
# - `name`: The name of the `Spring3` element (default: "").
# - `k`: The cubic spring stiffness constant (default: 1.0).
# - `x`: The initial displacement of the spring (default: 0.0).
# 
# ### Standard Elements
# 
# The entire system's elements, including mass, damper, linear spring, and the newly defined `Spring3`, are initialized:

# - The parameter values
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

# The forcing term can be set as an expression or even a Julia function. In this example, we will initialize a generic forcing term that will be defined later before solving the ODE.
@parameters F(t)
@named f = Se(F)

# ## Model Construction
# 
# Then, we build the system by simply passing the one-port elements to the 1-junction function as arguments.
@named mdl = Junction1([-1, m], [-1, d], [-1, s], [-1, s3], f)

# ## Analysing the model
# 
# We can visualize the bond-graph connections by generating a graph plot of the model.
generate_graph(mdl)

# The bond-graph equations can be extracted as:
equations(expand_connections(mdl))
latexify(equations(expand_connections(mdl)))

# The model's equations above were not simplified and represent the bond-graphs equations directly derived from the connections. Then, we simplify the DAE above and obtain the ODE. The simplified set of equations is required for using `solve` from `DifferentialEquations.jl`
@named sys = simplifysys(mdl);

# Print system states
states(sys)

# Print system parameters
parameters(sys)

# Print the system's simplified equations
equations(sys)

# By comparison, we can see that the equations above are the same as the canonical form of the harmonic oscillator presented at the beginning.

# We can generate LaTeX code from the equations with the following code.
latexify(equations(sys))

# ## System Simulations
# 
# ### Unforced System
# 
# Now, we are going to consider only the unforced case. Therefore, $F(t) = 0$. Thus, we update the force function on the system.
sys_unforced = substitute(sys, Dict(f.F => 0.0))

# Define the simulation time.
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

# When defining the ODEProblem as in the ModelingToolkit, we can define the initial value of the states.
prob = ODEProblem(sys_forced, [m.power.f => 0.0, s.q => 0.0], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")