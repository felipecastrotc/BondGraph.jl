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
@named m1 = Junction1([-1, m], [-1, d], [-1, s])

# nDoF spring and damper
@named sd = Junction1([-1, s], [-1, d])

# Set the number of mass
n = 4
# Generate the models for the masses
M, cons = [], []
for i = 2:n
    j0 = Junction0(sd; name = Symbol("m" * string(i) * "_j0"))
    j1 = Junction1([-1, m]; name = Symbol("m" * string(i) * "_j1"))
    push!(M, [j0, j1])
    push!(cons, connect(j0.power, j1.power))
end
# Connect the masses
push!(cons, [connect(M[i][2].power, M[i+1][1].power) for i = 1:(length(M)-1)]...)

# Connect first mass to the second
push!(cons, connect(m1.power, M[1][1].power))

# Flatten masses array
M = vcat(M...);

# Build the system
@named mdl = ODESystem(cons, t)
mdl = compose(mdl, m1, M...)

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
prob = ODEProblem(sys, [m1.s.q => 1.0], (0, 10))
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")
