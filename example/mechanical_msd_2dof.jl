using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify
using JLD2

# Define the elements
@named m = Mass(m = 1.0)      # (kg)
@named s = Spring(k = 1.0)    # (N/m)
@named d = Damper(c = 1.0)    # (N⋅m⋅s)

# Ground connected mass, mass 1
@named m1_j1 = Junction1([-1, m], [-1, d], [-1, s])

# 2nd DoF, mass 2
@named m2_j1 = Junction1([-1, m])
@named m2_j0 = Junction0()
@named m2_sd = Junction1([-1, s], [-1, d])

# Define the connections
cons = [
    connect(m2_j0.power, m2_j1.power),
    connect(m2_j0.power, m2_sd.power),
    connect(m1_j1.power, m2_j0.power),
]

# Build the system
@named mdl = ODESystem(cons, t)
mdl = compose(mdl, m1_j1, m2_j0, m2_j1, m2_sd)

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
prob = ODEProblem(sys, [m2_j1.m.power.f => 0.0, m2_sd.s.q => 1.0], (0, 10))
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")
