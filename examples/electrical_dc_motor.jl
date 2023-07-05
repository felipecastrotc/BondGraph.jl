# # [DC Motor](@id dc-motor)

#-
using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

# Reference
# https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling

# -----------------------------------------------------------------------------
# Implicit connection way

# Define the elements
@named L = Mass(m = 0.5)      # (H) Electric inductance
@named R = Damper(c = 1.0)    # (Ohm) Electric resistance
@named Uₐ = Se(12.0)        # (V) Voltage

@named J = Mass(m = 0.01)     # (kg⋅m^2) Shaft moment of inertia
@named b = Damper(c = 0.1)    # (N⋅m⋅s) Shaft viscous friction
@named Tₗ = Se(1.0)          # (N⋅m) Applied torque on the shaft

g = 0.01    # DC motor torque constant

# Mechanical part
@named jm = Junction1(Tₗ, [-1, b], [-1, J])
# Electical part
@named je = Junction1(Uₐ, [-1, R], [-1, L])

# Couple mechanical and electrical part
eqs = []
@named gy = mGY(je, jm, g = g, coneqs = eqs)

# Build the system
@named mdl = ODESystem(eqs, t)
mdl = compose(mdl, gy, je, jm)

generate_graph(mdl)
# Expand and simplify the system
@named sys = simplifysys(mdl)

# Print states, parameters, and equations
states(sys)
parameters(sys)
equations(sys)
# Generate LaTeX code from the equations
latexify(equations(sys))

tspan = (0, 10)
prob = ODEProblem(sys, [], tspan)
sol = solve(prob)
plot(sol, xlabel = "Time", ylabel = "Amplitude")

# -----------------------------------------------------------------------------
# Theoretical model

@variables i(t) = 0.0
@variables θ(t) = 0.0
@variables θ̇(t) = 0.0

Uᵥ = 12.0   # (V) - Voltage
Lᵥ = 0.5    # (H) - Inductance
Rᵥ = 1.0    # (Ohm) - Resistance
Jᵥ = 0.01   # (kg*m²) - Moment of intertia of the rotor
bᵥ = 0.1    # (N*m*s) - Motor viscous friction constant
Tᵥ = 1.0    # (N*m) - Motor torque

eqsₚ = [Jᵥ * D(θ̇) ~ Tᵥ + g * i - bᵥ * θ̇, Lᵥ * D(i) ~ -Rᵥ * i + Uᵥ - g * θ̇, D(θ) ~ θ̇]
@named sysₚ = ODESystem(eqsₚ, t)
sysₚ = structural_simplify(sysₚ)

probₚ = ODEProblem(sysₚ, [], tspan)
solₚ = solve(probₚ, reltol = 1e-8, abstol = 1e-8)

plot(solₚ.t, solₚ[θ̇])
plot!(solₚ.t, solₚ[i])


# -----------------------------------------------------------------------------
# Compare models

# Electrical part
plot(sol.t, sol[je.L.power.f], label = "Bond-graph", ylabel = "Current (A)")
plot!(solₚ.t, solₚ[i], label = "Model", xlabel = "Time (s)")

# Mechanical part
plot(sol.t, sol[jm.J.power.f], label = "Bond-graph", ylabel = "Angular velocity (rad/s)")
plot!(solₚ.t, solₚ[θ̇], label = "Model", xlabel = "Time (s)")
