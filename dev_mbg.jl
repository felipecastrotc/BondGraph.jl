using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")

# -----------------------------------------------------------------------------
# mTF - translation - CHAPTER 6. PLANAR MECHANICS -> master thesis zimmer_ms
@variables rx ry θ(t) = 0.0 a

d = [rx, ry]
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

@named m = Mass2D(m = [1.0, 1.0, 0.0])
nocon = ODESystem(Equation[], t, name = :nc)
# Sg = ODESystem2D(Se(0, name = :gx), Se(-9.81, name = :gy), Se(0, name = :gθ))

Sw = ODESystem2D(Sf(0, name = :wx), Sf(0, name = :wy), Se(0, name = :wθ))
Sg = ODESystem2D(nocon, Se(-9.81, name = :gy), nocon)
Sq = ODESystem2D(nocon, nocon, Dq(name = :qθ))

# @named j1 = Junction12D(Sw, -Sq)
@named j1 = Junction12D(Sw)
@named jm = Junction12D(Sg, -m)

# @named mtf = mTF2Dtrans(j1, jm, a = a, θ = θ, d = d)
θ = GlobalScope(θ)
# @named mtf = mTF2Dtrans(j1, jm, a = sqrt(2), θ = θ, d = [1.0, 1.0])
@named mtf = mTF2Dtrans(j1, jm, a = 1, θ = θ, d = [0.0, -1.0])

equations(mtf.x)
equations(mtf.y)
equations(mtf.θ)

@named mdl = get2Dsystem(mtf)

# eqs = [θ ~ Sq.θ.q]
eqs = [D(θ) ~ jm.θ.mθ.f]
mdl = extend(ODESystem(eqs, t, [], []; name = :mdl), mdl)

# sys = reducedobs(sys, name=:mdl)
sys = structural_simplify(mdl)

equations(sys)

prob = ODEProblem(sys, [], (0.0, 4.0))
sol = solve(prob)
plot(sol)

# equations(sys)
# equations(reducedobs(sys, name = :tst))

# sys = expand_connections(mdl)
# sys = alias_elimination(sys)
# sys2 = ode_order_lowering(sys)
# sys3 = dae_index_lowering(sys2)

# equations(sys)
# equations(sys2)
# equations(sys3)

# observed(sys2)
# equations(sys2)

# @named sys3 = reducedobs(sys2)

# equations(sys3)