using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D

using DifferentialEquations

# =============================================================================
# Function mGY dev place

function mGY(subsys...; name, g = 1.0)

    # Get connections
    c = subsys isa ODESystem ? [subsys, nothing] : collect(subsys)

    # Remove nothing from c array
    pos = .!isnothing.(c)
    c = c[pos]
    @assert !isempty(c)

    # If only one subsys is passed it automatically generates an open  
    # connection
    if sum(pos) == 1
        @named power = Power(type=tp)
        # Set variables according to the position
        if pos[1]
            e₁, f₁ = power.e, power.f
            e₂, f₂ = c[1].power.e, c[1].power.f
        else
            e₂, f₂ = power.e, power.f
            e₁, f₁ = c[1].power.e, c[1].power.f
        end
    else
        @assert length(c) == 2
        e₁, f₁ = c[1].power.e, c[1].power.f
        e₂, f₂ = c[2].power.e, c[2].power.f
    end

    # Gyrator equation
    eqs = [e₂ ~ g * f₁, e₁ ~ g * f₂]

    # Check if it is a modulated gyrator
    if isvariable(g)
        sts, ps = [], [g]
    elseif istree(unwrap(g))
        sts = []
        ps = collect(Set(ModelingToolkit.get_variables(g)))
    else
        sts, ps = [], []
    end

    sys = ODESystem(eqs, t, sts, ps; name = name)

    if @isdefined power
        return sys = compose(sys, power, c...)
    else
        return sys = compose(sys, c...)
    end
end

# =============================================================================
# DC motor
# -----------------------------------------------------------------------------
# Setting equations by hand

@named L = Mass(m = 0.5)
@named R = Damper(c = 1.0)
@named J = Mass(m = 0.01)
@named b = Damper(c = 0.1)

@variables e[1:8](t), f[1:8](t)
@parameters V, k, T

g = 0.01
g = k

sts = [f[4], f[5], e[1], e[4], e[5], e[8]]
eqs = [
    0 ~ e[1] - R.e - L.e - e[4],
    R.f ~ L.f,
    L.f ~ f[4],
    e[5] ~ g * f[4],
    e[4] ~ g * f[5],
    0 ~ e[5] - b.e - J.e - e[8],
    f[5] ~ b.f,
    b.f ~ J.f,
    # e[1] ~ V,
    # e[8] ~ T,
    e[1] ~ 12.0,
    e[8] ~ -1.0,
]

# g = 0.01
# g = k
mdl = compose(ODESystem(eqs, t, sts, []; name = :mdl), L, R, J, b)
equations(mdl)

sys = structural_simplify(mdl)
equations(sys)
states(sys)
parameters(sys)

prob = ODEProblem(sys, [], (0.0, 4.0))
sol = solve(prob)
plot(sol)

@named sysr = reducedobs(sys)

equations(sysr)
states(sysr)
parameters(sysr)

prob = ODEProblem(sysr, [], (0.0, 4.0))
sol = solve(prob)
plot(sol)

# -----------------------------------------------------------------------------
# Setting equations by functions

@named L = Mass(m = 0.5)
@named R = Damper(c = 1.0)
@named Uₐ = Se(12.0)

@named J = Mass(m = 0.01)
@named b = Damper(c = 0.1)
@named Tₗ = Se(1.0)

g = 0.01

@named je = Junction1(Uₐ, -R, -L, sgn = -1)
@named jm = Junction1(Tₗ, -b, -J)
@named gy = mGY(je, jm, g = g)

equations(gy)
@named sys = reducedobs(structural_simplify(gy))
equations(sys)

prob = ODEProblem(sys, [], (0.0, 4.0))
sol = solve(prob)
plot(sol)

# Example coherent with literature
# https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling

# -----------------------------------------------------------------------------
# Theory model

@variables i(t) = 0.0
@variables θ(t) = 0.0
@variables θ̇(t) = 0.0

Uᵥ = 12.0   # V - Voltage
Lᵥ = 0.5    # H - Inductance
Rᵥ = 1.0    # Ohm - Resistance
Jᵥ = 0.01   # kg*m² - Moment of intertia of the rotor
bᵥ = 0.1    # N*m*s - Motor viscous friction constant
Tᵥ = 1.0    # N*m - Motor torque

eqsₚ = [Jᵥ * D(θ̇) ~ Tᵥ + g * i - bᵥ * θ̇, Lᵥ * D(i) ~ -Rᵥ * i + Uᵥ - g * θ̇, D(θ) ~ θ̇]
@named sysₚ = ODESystem(eqsₚ, t)
sysₚ = structural_simplify(sysₚ)

probₚ = ODEProblem(sysₚ, [], (0.0, 4.0))
solₚ = solve(probₚ, reltol = 1e-8, abstol = 1e-8)
plot(solₚ.t, solₚ[θ̇])
plot(solₚ.t, solₚ[i])

# ------------------------------------------------------------------------------
# Comparison Solutions
solₚ = solve(probₚ, reltol = 1e-8, abstol = 1e-8)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

# Configure plot
ts = 0:0.01:4

stₚ = solₚ(ts)
st = sol(ts)

plot(stₚ.t, stₚ[2, :], xlabel = "Time (s)", label = "Theoretical")
plot!(st.t, st[2, :], line = (:dot, 4), ylabel = "I (A)", label = "BG mGY")
plot!(size = 72 .* (5.5, 3.5), dpi = 300)
savefig("/home/fctc/vm-share/Equinor/images/dc_sim_A.png")


plot(stₚ.t, stₚ[3, :], xlabel = "Time (s)", label = "Theoretical")
plot!(st.t, st[1, :], line = (:dot, 4), ylabel = "θ̇ (rad/s)", label = "BG mGY")
plot!(size = 72 .* (5.5, 3.5), dpi = 300)
savefig("/home/fctc/vm-share/Equinor/images/dc_sim_rad.png")