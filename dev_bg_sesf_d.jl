using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")

# =============================================================================
# Function Se

function Se(expr; name)
    sts = @variables e(t) = 0.0
    eqs = [sts[1] ~ expr]
    ODESystem(eqs, t, sts, []; name = name)
end

function Sf(expr; name)
    sts = @variables f(t) = 0.0
    eqs = [sts[1] ~ expr]
    ODESystem(eqs, t, sts, []; name = name)
end

function Dq(; name, x = 0.0)
    @variables f(t) = 0.0
    @variables q(t) = 0.0

    eqs = [
        D(q) ~ f
    ]
    ODESystem(eqs, t, [q, f], []; name = name)
end

# -----------------------------------------------------------------------------
# Dev

@variables θ
θ = GlobalScope(θ)

@named sf = Sf(sin(θ))
@named se = Se(sin(θ))
@named dq = Dq()
dq.q

# Test junciton1
@named m = Mass(m = 0.5)
@named s = Spring(k = 1.0)
@named d = Damper(c = 1.0)

# Simple MSD system
@named c = Junction1(m, s, d, se, dq, couple = false)

# Add the θ equation
eqs = [θ ~ dq.q]
c = extend(ODESystem(eqs, t, [], []; name = :c), c)

# Simplify
equations(c)

sys = structural_simplify(c)
@named sys = reducedobs(sys)
equations(sys)