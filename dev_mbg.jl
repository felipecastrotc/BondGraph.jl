using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")

# -----------------------------------------------------------------------------
# mTF - translation - CHAPTER 6. PLANAR MECHANICS -> master thesis zimmer_ms

@variables rx ry θ a

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
Sw = ODESystem2D(Sf(0, name = :wx), Sf(0, name = :wy), Se(0, name = :wθ))
Sg = ODESystem2D(Se(0, name = :gx), Se(-9.81, name = :gy), Se(0, name = :gθ))

@named j1 = Junction12D(Sw)
@named jm = Junction12D(Sg, -m)

# @named mtf = mTF2Dtrans(j1, jm, a = a, θ = θ, d = d)
@named mtf = mTF2Dtrans(j1, jm, a = 1, θ = θ, d = [0.0, 1.0])

equations(mtf.x)
equations(mtf.y)
equations(mtf.θ)

@named mdl = get2Dsystem(mtf)

eqs = [θ ~ j1.θ.f]
mdl = extend(ODESystem(eqs, t, [], []; name = :mdl), mdl)

structural_simplify(mdl)