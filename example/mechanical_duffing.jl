using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# https://en.wikipedia.org/wiki/Duffing_equation
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

# Set the parameters values
α = 1.0
β = 5.0
δ = 0.02
γ = 8.0
ω = 0.5
@parameters f(t)

# Define the elements
@named m = Mass(m = 1.0)
@named s = Spring(k = α)
@named s3 = Spring3(k = β)
@named d = Damper(c = δ)

# Force function
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
