using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D, mGY, tpgy
import ModelingToolkit: isvariable, istree
import Symbolics: unwrap, wrap

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

    # Generate the in and out connection
    @named pin = Power(type=tpgy)
    @named pout = Power(type=tpgy)

    # Alias for simpler view about the mGY port
    e₁, f₁ = pin.e, pin.f
    e₂, f₂ = pout.e, pout.f

    # Gyrator equation
    eqs = [e₂ ~ g * f₁, e₁ ~ g * -f₂]
    # TODO: check why adding the - on f₂ solves the signal issue on DC motor
    # TODO: I already tried to change the signal on the port type tpgy at 
    #       adjmtx2eqs

    # Check if it is a modulated gyrator
    if isvariable(g)
        sts, ps = [], [g]
    elseif istree(unwrap(g))
        vars = collect(Set(ModelingToolkit.get_variables(g)))
        sts = filter(x -> ~isindependent(Num(x)), vars)
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        sts, ps = [], []
    end

    # Apply connections
    if sum(pos) == 1
        if pos[1]
            push!(eqs, connect(c[1].power, pin))
        else
            push!(eqs, connect(pout, c[1].power))
        end
    else
        @assert length(c) == 2
        push!(eqs, connect(c[1].power, pin), connect(pout, c[2].power))
    end

    return compose(ODESystem(eqs, t, sts, ps; name = name), pin, pout)
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
@named Uₐ = Se(-12.0)

@named J = Mass(m = 0.01)
@named b = Damper(c = 0.1)
@named Tₗ = Se(-1.0)

g = 0.01

@named jm = Junction1(Tₗ, b, J)
@named je = Junction1(Uₐ, R, L)
# @named jm = Junction1(b, J)
@named gy = mGY(je,jm, g = g)
equations(gy)

@named mdl = compose(gy, je, jm)
equations(mdl)
generate_graph(mdl)
emdl = expand_connections(mdl)
@named sys = reducedobs(structural_simplify(emdl))
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
plot!(solₚ.t, solₚ[i])

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