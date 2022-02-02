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

# function Sf(expr; name)
#     @named power = Power(effort = nothing, flow = expr)
#     extend(ODESystem(Equation[], t, [], []; name = name), power)
# end

# -----------------------------------------------------------------------------
# Dev

@variables θ
θ = GlobalScope(θ)

@named sf = Sf(sin(θ))
@named se = Se(sin(θ))

# Test junciton1
@named m = Mass(m = 0.5)
@named s = Spring(k = 1.0)
@named d = Damper(c = 1.0)

# Simple MSD system
@named c = Junction1(m, s, d, se, couple = false)

# Add the θ equation
eqs = [θ ~ s.q]
c = extend(ODESystem(eqs, t, [], []; name = :c), c)

# Simplify
sys = structural_simplify(c)
@named sys = reducedobs(sys)
equations(sys)