using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# Define the elements
@named m = Mass(m = 1.0)      # (kg)
@named s = Spring(k = 1.0)    # (N/m)
@named d = Damper(c = 1.0)    # (N⋅m⋅s)

# Simple connection using the Junction
@named mdl = Junction1(m, d, s)

# Generate visual graph
generate_graph(mdl)

# Expand and simplify the system
@named sys = simplifysys(mdl)

# Print states, parameters, and equations
states(sys)
parameters(sys)
equations(sys)
# Generate LaTeX code from the equations
latexify(equations(sys))

# Solve the ODE problem
prob = ODEProblem(sys, [m.power.f => 0.0, s.q => 1.0], (0, 10))
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")
