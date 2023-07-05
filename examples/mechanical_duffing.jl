# # [Duffing Equation](@id msd-duffing)
# 
# The Duffing equation is an extension of the simple mass-spring-damper system that incorporates a cubic spring. Despite its simple form, the Duffing equation exhibits chaotic behavior, and it is used to model more realistic damped oscillators. The equation is given by:
# 
# $$m\,\frac{d^2 x}{dt^2} + c\,\frac{d x}{dt} + k\,x + k_3\,x^3 = F(t),$$
# 
# where $x$ represents the position, $m$ is the mass, $c$ is the damping coefficient, $k$ is the spring stiffness, $k_3$ is the cubic spring stiffness, and $F(t)$ denotes an external force applied to the system.
# 
# In standard bond graphs, the capacitance element is linear, while for the Duffing equation we need a non-linear capacitance, cubic spring. Therefore, to model the Duffing equation using the Bond Graph Toolkit, we need to introduce a non-linear capacitance element. The Bond Graph Toolkit allows users to define custom elements that utilize the `Power()` connector.
# 
# ## Bond Graph Toolkit
# 
# To begin, we import the Bond Graph Toolkit module for modeling the system and the `DifferentialEquations.jl` package for solving the resulting ordinary differential equation (ODE). Additionally, we import the independent variable `t` from the Bond Graph Toolkit, which will be used to define a custom forcing term.

using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# ## Building the model
# 
# ### Custom capacitance

# Duffing equation example

# Define the the cubic spring
function Spring3(; name, k = 1.0, x = 0.0)
    @named power = Power()

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        power.e ~ q^3 / C
        D(q) ~ power.f
    ]
    compose(ODESystem(eqs, t, [q], ps; name = name), power)
end

# ### Setting the elements
# 
# Prior to the definition of a system, we need to define the single port elements of the system. Then, we define the mass, spring and damper elements as follows:

# Set the parameters values
α = 1.0
β = 5.0
δ = 0.02
γ = 8.0
ω = 0.5

# Define the elements
@named m = Mass(m = 1.0)
@named s = Spring(k = α)
@named s3 = Spring3(k = β)
@named d = Damper(c = δ)

# Force function
@parameters f(t)
@named F = Se(f)

# Simple connection using the Junction
@named mdl = Junction1(m, d, s, s3, [-1, F])

# Generate visual graph
generate_graph(mdl)

# Expand and simplify the system
@named sys = simplifysys(mdl)

# Set the force function
sys = renamevars(sys, Dict(F.f => γ * cos(ω * t)))

# Print states, parameters, and equations
states(sys)
parameters(sys)
equations(sys)
# Generate LaTeX code from the equations
latexify(equations(sys))

# Solve the ODE problem
tspan = (0, 10)
prob = ODEProblem(sys, [], tspan)
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")
