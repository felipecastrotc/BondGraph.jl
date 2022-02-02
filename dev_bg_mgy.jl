using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")

# =============================================================================
# Function mGY

function mGY(subsys...; name, g = 1.0, couple = true)
    # Get connections
    c = collect(Base.Flatten([subsys]))

    if couple
        pos = .!isnothing.(c)
        @assert length(sum(pos)) == 1

        @named power = Power()

        # Set variables according to the position
        if pos[1]
            @unpack e₁, f₁ = power
            e₂, f₂ = c[1].e, c[1].f
        else
            @unpack e₂, f₂ = power
            e₁, f₁ = c[1].e, c[1].f
        end

    else
        @assert length(c) == 2

        e₁, f₁ = c[1].e, c[1].f
        e₂, f₂ = c[2].e, c[2].f
    end

    # Gyrator equation
    eqs = [
        e₂ ~ g * f₁,
        e₁ ~ g * f₂,
    ]

    # Remove nothing from c array
    filter!(!isnothing, c)

    # Check if it is a modulated gyrator
    if isvariable(g)
        sts, ps = [], [g]
    elseif istree(unwrap(g))
        sts = []
        ps = collect(Set(extract_vars(g)))
    else
        sts, ps = []
    end

    sys = ODESystem(eqs, t, sts, ps; name = name)

    if @isdefined power
        sys = extend(sys, power)
    end

    compose(sys, c...)
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
    e[8] ~ -1.,
]


# mdl = compose(ODESystem(eqs, t, sts, [V, k, T]; name = :mdl), L, R, J, b)
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
@named J = Mass(m = 0.01)
@named b = Damper(c = 0.1)

@parameters V, k, T






# mdl = compose(ODESystem(eqs, t, sts, [V, k, T]; name = :mdl), L, R, J, b)
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

