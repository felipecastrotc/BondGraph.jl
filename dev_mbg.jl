using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")

# -----------------------------------------------------------------------------
# mTF - translation - CHAPTER 6. PLANAR MECHANICS -> master thesis zimmer_ms
@variables rx ry θ(t) = 0.0 a θ₀

d = [rx, ry]
d = [1, 0]
A = [cos(θ) sin(θ); -sin(θ) cos(θ)] * d

[1 0 A[2]; 0 1 -A[1]; 0 0 1]

p = Dict("a" => a, "θ" => θ, "d" => d)
M = M2Dtrans(p)

# -----------------------------------------------------------------------------
# Test junctions

@named m = Mass2D()
@named s = Spring2D()
@named d = Damper2D()

@named j1 = Junction12D(m, s, d, couple = false)
@named j0 = Junction12D(m, s, d, couple = false)

# -----------------------------------------------------------------------------
# Test Pendulum Model

@named m = Mass2D(m = [2, 2, 2])
nocon = ODESystem(Equation[], t, name = :nc)

Sw = ODESystem2D(Sf(0, name = :wx), Sf(0, name = :wy), Se(0, name = :wθ))
Sg = ODESystem2D(Se(0, name = :gx), Se(-9.81, name = :gy), Se(0, name = :gθ))
# Sg = ODESystem2D(nocon, Se(-9.81*2, name = :gy), nocon)
Sq = ODESystem2D(nocon, nocon, Dq(name = :qθ))

# @named j1 = Junction12D(Sw, -Sq)
@named j1 = Junction12D(Sw)
@named jm = Junction12D(Sg, -m)

θ = GlobalScope(θ)
A = 1.0
θ₀ = deg2rad(0.1)
Rx = A * sin(θ₀)
Ry = A * cos(θ₀)
@named mtf = mTF2Dtrans(j1, jm, a = A, θ = θ, d = [Rx, Ry])
@named mdl = get2Dsystem(mtf)

eqs = [D(θ) ~ jm.θ.mθ.f]
mdl = extend(ODESystem(eqs, t, [], []; name = :mdl), mdl)

sys = structural_simplify(mdl)
sys = reducedobs(sys, name=:mdl)
equations(sys)

prob = ODEProblem(sys, [], (0.0, 20.0))
sol = solve(prob)
# plot(sol.t, sqrt.(sol[jm.y.e] .^ 2 + sol[jm.x.e] .^ 2))
plot(sol)
plot(sol.t, sol[θ])

