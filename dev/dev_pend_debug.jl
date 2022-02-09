using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D

using DifferentialEquations


# Global definitions
θ₀ = 0.999 * π
l₀ = 10
m₀ = 1

@parameters l, g
@variables θ(t) = θ₀
g = GlobalScope(g)
l = GlobalScope(l)
θ = GlobalScope(θ)

# -----------------------------------------------------------------------------
# Theory model
@variables θ̇(t) = 0.0

eqsₚ = [D(θ̇) ~ -(g / l) * sin(θ), D(θ) ~ θ̇]
@named sysₚ = ODESystem(eqsₚ, t)
sysₚ = structural_simplify(sysₚ)
equations(sysₚ)

probₚ = ODEProblem(sysₚ, [θ => θ₀], (0.0, 30.0), [l => l₀, g => 9.81])
sol = solve(probₚ, reltol = 1e-8, abstol = 1e-8)
plot(sol)

# -----------------------------------------------------------------------------
# Bond 1D non-linear compliance

function SpringSin(; name, m = 1.0, g = 9.81, l = 1.0, θ = 0.0)
    @named power = Power()
    @unpack e, f = power

    @variables q(t) = θ
    ps = @parameters g = g l = l

    eqs = [
        e ~ sin(q) * g * m / l,
        D(q) ~ f
    ]
    extend(ODESystem(eqs, t, [q], ps; name = name), power)
end

@named m₁ = Mass(m = m₀)
@named s₁ = SpringSin(m = m₀, g = g, l = l, θ = θ)
@named p₁ = Junction1(m₁, s₁, couple = false)

@named sys₁ = reducedobs(structural_simplify(p₁))
equations(sys₁)

prob₁ = ODEProblem(sys₁, [θ => θ₀], (0.0, 40.0), [l => l₀, g => 9.81])
sol₁ = solve(prob₁, reltol = 1e-8, abstol = 1e-8)
plot(sol₁)


# -----------------------------------------------------------------------------
# Yellow blue Using Junctions and mTR

@named mⱼ = Mass(m = m₀ * (1 - 0.0001))
@named Jⱼ = Mass(m = m₀ * 0.0001)
@named gⱼ = Se(9.81 * m₀)

@named x = Junction1(-mⱼ)
@named xtf = mTF(x, r = (l₀ * cos(θ)))

@named y = Junction1(-mⱼ, gⱼ)
@named ytf = mTF(y, r = -(l₀ * sin(θ)))
@named mdlⱼ = Junction1(-Jⱼ, -ytf, -xtf, couple = false)

eqsⱼ = [D(θ) ~ Jⱼ.f]
mdl2ⱼ = extend(ODESystem(eqsⱼ, t, [θ], []; name = :mdl), mdlⱼ)
@named sysⱼ = reducedobs(structural_simplify(mdl2ⱼ))
equations(structural_simplify(mdl2ⱼ))

probⱼ = ODEProblem(sysⱼ, [θ => θ₀], (0.0, 40.0))
solⱼ = solve(probⱼ, reltol = 1e-8, abstol = 1e-8)
# plot(solⱼ)
# plot(solⱼ.t, solⱼ[θ])
# plot!(solⱼ.t, solⱼ[Jⱼ.f])

plot(sol₁.t, sol₁[2, :], label = "target")
plot!(solⱼ.t, solⱼ[θ], label = "sim")
# ------------------------------------------------------------------------------
# Comparison Solutions
solₚ = solve(probₚ, reltol = 1e-8, abstol = 1e-8)
sol₁ = solve(prob₁, reltol = 1e-8, abstol = 1e-8)
solⱼ = solve(probⱼ, reltol = 1e-8, abstol = 1e-8)

# Configure plot
ts = 0:0.01:25

stₚ = solₚ(ts)
st₁ = sol₁(ts)
stⱼ = solⱼ(ts)
plot(stₚ.t, stₚ[2, :], xlabel = "Time (s)", label = "Theoretical")
plot!(st₁.t, st₁[2, :], ylabel = "θ (rad)", label = "BG Custom")
plot!(stⱼ.t, stⱼ[θ], label = "BG mTF")

plot(stₚ.t, stₚ[1, :], xlabel = "Time (s)", label = "Theoretical")
plot!(st₁.t, st₁[1, :], ylabel = "θ (rad/s)", label = "BG Custom")
plot!(stⱼ.t, stⱼ[Jⱼ.f], label = "BG mTF")



Θ₀ = (π * 2 .^ (-3:0.2:0) .- 1e-4)
plot(title="Phase space", xlabel="θ", ylabel="θ vel", legend=false)
for i in Θ₀
    probₚ = ODEProblem(sysₚ, [θ => i], (0.0, 60.0), [l => l₀, g => 9.81])
    ϕ = solve(probₚ, reltol = 1e-8, abstol = 1e-8)
    plot!(ϕ[2, :], ϕ[1, :], lc=:blue)
end
plot!()

# -----------------------------------------------------------------------------
# Bond 2D

@parameters rx ry
@variables θ(t) = 0

θ = GlobalScope(θ)
rx, ry = GlobalScope(rx), GlobalScope(ry)

Sw = ODESystem2D(Sf(0, name = :wx), Sf(0, name = :wy), Se(0, name = :wθ))
Sg = ODESystem2D(Se(0, name = :gx), Se(-9.81, name = :gy), Se(0, name = :gθ))

@named m = Mass2D(m = [1, 1, 0.0001])

@named j1 = Junction12D(Sw)
@named jm = Junction12D(Sg, -m)

A = l₀
@named mtf = mTF2Dtrans(j1, jm, a = A, θ = θ, d = [rx, ry])
@named mdl = get2Dsystem(mtf)

eqs = [D(θ) ~ jm.θ.mθ.f]
mdl = extend(ODESystem(eqs, t, [], []; name = :mdl), mdl)
sys = structural_simplify(mdl)

# θ₀ = deg2rad(0.1)
θ₀₂ = π - θ₀
p0 = [θ => pi, rx => A * sin(θ₀₂), ry => A * cos(θ₀₂)]
prob₂ = ODEProblem(sys, [], (0.0, 80.0), p0)
sol = solve(prob₂, reltol = 1e-8, abstol = 1e-8)
plot(sol)
plot(sol.t, sol[θ])

# Backcup
st₂ = sol₂(ts)
plot!(st₂.t, -st₂[θ] .+ pi)